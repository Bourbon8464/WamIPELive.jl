module WamIPELive

using Dates
using Printf
using HTTP
using URIs
using EzXML
using NCDatasets
using Serialization
using Statistics
using DataInterpolations

export WFSInterpolator, get_value, get_batch, print_cache_stats, list_vars, dump_sample, dump_all

# =============================== Configuration ===============================

"""
    WFSInterpolator; kwargs...

Configuration for the NOMADS WAM–IPE v1.2 HTTP client and spatio-temporal interpolation.

# Keyword Arguments
- `base_url::String="https://nomads.ncep.noaa.gov/pub/data/nccf/com/wfs/v1.2"`
- `product::String="wfs"`: `"wfs"` or `"wrs"`.
- `stream::String="ipe10"`: `"ipe05"`, `"ipe10"`, `"gsm05"`, `"gsm10"`.
- `varname::String="ion_temperature"`: NetCDF variable to sample.
- `interpolation::Symbol=:sciml`: `:nearest`, `:linear`, `:logz_linear`, `:logz_quadratic`, or `:sciml` (alias for quadratic in `log(z)` with bilinear lon/lat).
- `cache_dir::String=normpath("./cache")`
- `cache_max_bytes::Int=2_000_000_000`
"""
Base.@kwdef struct WFSInterpolator
    base_url::String = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/wfs/v1.2"
    product::String = "wfs"
    stream::String = "ipe10"
    varname::String = "ion_temperature"
    interpolation::Symbol = :sciml
    cache_dir::String = normpath("./cache")
    cache_max_bytes::Int = 2_000_000_000
end

_normalize_interp(s::Symbol) = (s === :sciml ? :logz_quadratic : s)
const _ALLOWED_INTERP_NORM = Set([:nearest, :linear, :logz_linear, :logz_quadratic])
const _UA = "WFS.jl/0.1 (+https://example.invalid)"
const _HDRS = ["User-Agent" => _UA, "Accept" => "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8"]

# ================================ HTTP layer =================================

function _http_get_classify(url::AbstractString; readtimeout::Int=60)
    try
        r = HTTP.get(url; readtimeout=readtimeout, headers=_HDRS)
        return r.status == 200 ? r : r.status == 404 ? :notfound : r.status == 403 ? :forbidden :
               error("HTTP $(r.status) for $url")
    catch e
        if e isa HTTP.Exceptions.StatusError
            s = e.status
            return s == 404 ? :notfound : s == 403 ? :forbidden : rethrow()
        else
            rethrow()
        end
    end
end

function _http_head_exists(url::AbstractString; readtimeout::Int=30)::Union{Bool,Symbol}
    try
        r = HTTP.request("HEAD", url; headers=_HDRS, readtimeout=readtimeout)
        return r.status == 200 ? true : r.status == 404 ? false : r.status == 403 ? :forbidden : false
    catch
        try
            r = HTTP.request("GET", url; headers=vcat(_HDRS, ["Range"=>"bytes=0-0"]), readtimeout=readtimeout)
            return r.status in (200,206) ? true : r.status == 404 ? false : r.status == 403 ? :forbidden : false
        catch
            return false
        end
    end
end

# ================================== Cache ====================================

const _CACHE_META_FILE = "metadata.bin"
mutable struct _FileCache
    dir::String
    max_bytes::Int64
    map::Dict{String,String}
    sizes::Dict{String,Int64}
    order::Vector{String}
    bytes::Int64
    lock::ReentrantLock
end
const _CACHES = Dict{Tuple{String,Int64}, _FileCache}()

function _load_cache(dir::AbstractString, max_bytes::Int64)
    mkpath(dir)
    meta = joinpath(dir, _CACHE_META_FILE)
    if isfile(meta)
        try
            open(meta,"r") do io
                c = deserialize(io)
                if c isa _FileCache
                    c.bytes = sum(values(c.sizes))
                    c.order = [k for k in c.order if haskey(c.map,k) && isfile(c.map[k])]
                    c.lock = ReentrantLock()
                    return c
                end
            end
        catch
        end
    end
    _FileCache(String(dir), max_bytes, Dict(), Dict(), String[], 0, ReentrantLock())
end
function _save_cache(c::_FileCache)
    open(joinpath(c.dir,_CACHE_META_FILE),"w") do io; serialize(io,c); end
end
function _get_cache(dir::AbstractString, max_bytes::Int64)
    key = (String(dir),Int64(max_bytes))
    get!(_CACHES, key) do; _load_cache(dir,max_bytes) end
end
function print_cache_stats(itp::WFSInterpolator)
    c = _get_cache(itp.cache_dir, itp.cache_max_bytes)
    lock(c.lock) do
        println("Cache dir: ",c.dir)
        println("Capacity : ",round(c.max_bytes/1e9,digits=2)," GB")
        println("Used     : ",round(c.bytes/1e9,digits=3)," GB (",length(c.map)," files)")
        if !isempty(c.order)
            println("LRU head: ",first(c.order))
            println("MRU tail: ",last(c.order))
        end
    end
end
function _touch!(c::_FileCache, key::String)
    idx = findfirst(==(key), c.order)
    idx !== nothing && deleteat!(c.order, idx)
    push!(c.order, key)
end
function _evict!(c::_FileCache)
    while c.bytes > c.max_bytes && !isempty(c.order)
        k = first(c.order); popfirst!(c.order)
        if haskey(c.map,k)
            path = c.map[k]; sz = get(c.sizes,k,0)
            try; isfile(path) && rm(path; force=true); catch; end
            delete!(c.map,k); delete!(c.sizes,k)
            c.bytes = max(0, c.bytes - sz)
        end
    end
end

function _download_http_cached(itp::WFSInterpolator, url::AbstractString)
    c   = _get_cache(itp.cache_dir, itp.cache_max_bytes)
    key = String(url)
    u   = URIs.URI(url)
    rel = lstrip(u.path, '/')
    local_path = normpath(joinpath(itp.cache_dir, rel))
    lock(c.lock) do
        if haskey(c.map, key) && isfile(c.map[key])
            _touch!(c, key); _save_cache(c)
            return c.map[key]
        end
    end
    mkpath(dirname(local_path))
    tmp = local_path * ".part"
    r = try
        HTTP.request("GET", url; headers=_HDRS, readtimeout=120)
    catch e
        error("GET failed for $url: $(e)")
    end
    if r.status == 200
        open(tmp, "w") do f; write(f, r.body); end
        mv(tmp, local_path; force=true)
        sz = filesize(local_path)
        lock(c.lock) do
            c.map[key]   = local_path
            c.sizes[key] = sz
            c.bytes     += sz
            _touch!(c, key); _evict!(c); _save_cache(c)
        end
        return local_path
    elseif r.status in (403,404)
        isfile(tmp) && try rm(tmp; force=true) catch; end
        error("HTTP $(r.status) for $url (blocked or missing)")
    else
        isfile(tmp) && try rm(tmp; force=true) catch; end
        error("HTTP $(r.status) for $url")
    end
end

# ============================ NOMADS directory ===============================

function _cycle_for_wfs(dt::DateTime)::DateTime
    h = hour(dt)
    if h < 3; DateTime(Date(dt), Time(0))
    elseif h < 9; DateTime(Date(dt), Time(6))
    elseif h < 15; DateTime(Date(dt), Time(12))
    elseif h < 21; DateTime(Date(dt), Time(18))
    else; DateTime(Date(dt)+Day(1), Time(0))
    end
end
function _cycle_for_wrs(dt::DateTime)::DateTime
    h = hour(dt)
    if h < 3;  DateTime(Date(dt)-Day(1), Time(18))
    elseif h < 9; DateTime(Date(dt), Time(0))
    elseif h < 15; DateTime(Date(dt), Time(6))
    elseif h < 21; DateTime(Date(dt), Time(12))
    else; DateTime(Date(dt), Time(18))
    end
end
function _cycle(itp::WFSInterpolator, dt::DateTime)
    itp.product == "wfs" ? _cycle_for_wfs(dt) :
    itp.product == "wrs" ? _cycle_for_wrs(dt) :
    error("product must be wfs or wrs")
end
function _dir_url_cycle(itp::WFSInterpolator, dt::DateTime)
    cyc = _cycle(itp, dt)
    ymd = Dates.format(Date(cyc), dateformat"yyyymmdd")
    HH  = @sprintf("%02d", hour(cyc))
    string(itp.base_url, "/", itp.product, ".", ymd, "/", HH, "/")
end
_parse_vtime_from_name(name::AbstractString) = let m = match(r"(\d{8})_(\d{6})\.nc$", name)
    m === nothing && return nothing
    ymd, hms = m.captures
    DateTime(parse(Int, ymd[1:4]), parse(Int, ymd[5:6]), parse(Int, ymd[7:8]),
             parse(Int, hms[1:2]), parse(Int, hms[3:4]), parse(Int, hms[5:6]))
end
function _dir_url_for_cycle_datetime(itp::WFSInterpolator, cyc::DateTime)
    ymd = Dates.format(Date(cyc), dateformat"yyyymmdd")
    HH  = @sprintf("%02d", hour(cyc))
    string(itp.base_url, "/", itp.product, ".", ymd, "/", HH, "/")
end

function _list_hrefs_safe(url::AbstractString)::Vector{String}
    r = _http_get_classify(url; readtimeout=60)
    r === :notfound  && return String[]
    r === :forbidden && return String[]
    r isa HTTP.Messages.Response || return String[]
    doc  = EzXML.parsehtml(String(r.body))
    hrefs = String[]
    for a in EzXML.findall("//a", doc.root)
        if haskey(a, "href")
            h = a["href"]
            if h !== nothing && endswith(h, ".nc")
                push!(hrefs, String(h))
            end
        end
    end
    return hrefs
end

_cadence_minutes(stream::AbstractString) = Base.startswith(lowercase(stream), "gsm05") ? 5 :
                                           Base.startswith(lowercase(stream), "gsm10") ? 10 :
                                           Base.startswith(lowercase(stream), "ipe05") ? 5 :
                                           Base.startswith(lowercase(stream), "ipe10") ? 10 : 10

function _build_filename(product::String, cyc_hour::Int, token::String, vt::DateTime)
    @sprintf("%s.t%02dz.%s.%04d%02d%02d_%02d%02d%02d.nc",
             product, cyc_hour, token,
             year(vt), month(vt), day(vt), hour(vt), minute(vt), second(vt))
end

function _normalize_stream(s::AbstractString)
    ls = lowercase(s)
    startswith(ls, "gsm") ? replace(ls, "gsm" => "wam") : ls
end

function _discover_files_by_probe(itp::WFSInterpolator, dt::DateTime; window_minutes::Int=180)
    cyc = _cycle(itp, dt)
    ymd = Dates.format(Date(cyc), dateformat"yyyymmdd")
    HH  = @sprintf("%02d", hour(cyc))
    base = string(itp.base_url, "/", itp.product, ".", ymd, "/", HH, "/")
    token = _normalize_stream(itp.stream)
    cad   = _cadence_minutes(itp.stream)
    base_vt = DateTime(floor(DateTime(dt), Minute(cad)))
    vs = DateTime[]
    for m in -window_minutes:cad:window_minutes
        push!(vs, base_vt + Minute(m))
    end
    found = Dict{DateTime,String}()
    for vt in vs
        fname = _build_filename(itp.product, hour(cyc), token, vt)
        url = base * fname
        ex = _http_head_exists(url)
        ex === true && (found[vt] = url)
    end
    return found
end

function _pick_two_files(itp::WFSInterpolator, dt::DateTime)
    base_cyc = _cycle(itp, dt)
    tok = _normalize_stream(itp.stream)
    cad = _cadence_minutes(itp.stream)
    for back in 0:8
        cyc = base_cyc - Hour(6*back)
        url_dir = _dir_url_for_cycle_datetime(itp, cyc)
        names   = _list_hrefs_safe(url_dir)
        if isempty(names)
            function probe_in_cycle(dt_probe::DateTime)
                found = Dict{DateTime,String}()
                ymd = Dates.format(Date(cyc), dateformat"yyyymmdd")
                HH  = @sprintf("%02d", hour(cyc))
                base = string(itp.base_url, "/", itp.product, ".", ymd, "/", HH, "/")
                base_vt = DateTime(floor(DateTime(dt_probe), Minute(cad)))
                for m in -180:cad:180
                    vt = base_vt + Minute(m)
                    fname = _build_filename(itp.product, hour(cyc), tok, vt)
                    url   = base * fname
                    ex = _http_head_exists(url)
                    ex === true && (found[vt] = url)
                end
                return found
            end
            found = probe_in_cycle(dt)
            if isempty(found)
                found = probe_in_cycle(DateTime(Date(cyc), Time(hour(cyc))))
            end
            if !isempty(found)
                vtimes = sort!(collect(keys(found)))
                ilo = searchsortedlast(vtimes, dt)
                ihi = searchsortedfirst(vtimes, dt)
                nlo = ilo == 0             ? vtimes[1]   : vtimes[ilo]
                nhi = ihi > length(vtimes) ? vtimes[end] : vtimes[ihi]
                return found[nlo], found[nhi], nlo, nhi
            else
                continue
            end
        end
        cand = [n for n in names if occursin(tok, lowercase(n))]
        isempty(cand) && (cand = names)
        times = Dict{String,DateTime}()
        for n in cand
            vt = _parse_vtime_from_name(n)
            vt === nothing && continue
            times[n] = vt
        end
        isempty(times) && continue
        vtq = dt
        before = [n for n in keys(times) if times[n] <= vtq]
        after  = [n for n in keys(times) if times[n] >= vtq]
        nlo = !isempty(before) ? before[argmax([times[n] for n in before])] :
                                 begin
                                     all = collect(keys(times))
                                     all[argmin([abs(times[n]-vtq) for n in all])]
                                 end
        nhi = !isempty(after)  ? after[argmin([times[n] for n in after])]  :
                                 begin
                                     all = collect(keys(times))
                                     all[argmin([abs(times[n]-vtq) for n in all])]
                                 end
        return string(url_dir, nlo), string(url_dir, nhi), times[nlo], times[nhi]
    end
    error("no files found by listing or probing recent cycles for $(itp.product) around $(dt)")
end

# ================================== Grids ====================================

_grid_uses_360(lon::AbstractVector) = maximum(lon) > 180
function _wrap_lon_for_grid(lon_grid::AbstractVector, lonq::Real)
    _grid_uses_360(lon_grid) ? (lonq < 0 ? lonq + 360 : lonq) : (lonq > 180 ? lonq - 360 : lonq)
end

function _classify_vertical_units(units_raw::AbstractString)
    s = lowercase(strip(String(units_raw)))
    isempty(s) && return :missing
    occursin(r"\bkm\b", s) && return :km
    if (occursin(r"\bm\b", s) || occursin("meter", s)) && !occursin("km", s); return :m; end
    (occursin("pa",s) || occursin("mb",s) || occursin("hpa",s) || occursin("pressure",s)) && return :pressure
    (occursin("level",s) || occursin("index",s) || occursin("layer",s)) && return :index
    :unknown
end

function _load_grids(ds::NCDataset, varname::String; file_time::Union{DateTime,Nothing}=nothing)
    haskey(ds, varname) || error("variable '$varname' not found")
    v        = ds[varname]
    dnamesT  = NCDatasets.dimnames(v)
    dnames   = collect(String.(dnamesT))
    Vraw     = Array(v)
    nd       = ndims(Vraw)
    function classify(d::String)
        ld    = lowercase(d)
        a     = haskey(ds, d) ? Dict(ds[d].attrib) : Dict{String,Any}()
        std   = lowercase(string(get(a, "standard_name", "")))
        axis  = uppercase(string(get(a, "axis", "")))
        units = lowercase(string(get(a, "units", "")))
        if occursin("time", ld) || axis == "T" || std == "time"; return :time end
        if occursin("lat",  ld) || std == "latitude"  || axis == "Y" || occursin("degrees_north", units); return :lat end
        if occursin("lon",  ld) || std == "longitude" || axis == "X" || occursin("degrees_east",  units); return :lon end
        if occursin("lev",  ld) || occursin("height", ld) || occursin("alt", ld) || ld == "z" || axis == "Z"; return :z end
        if ld in ("x","grid_xt","i","nx"); return :lon end
        if ld in ("y","grid_yt","j","ny"); return :lat end
        :unknown
    end
    roles = collect(map(classify, dnames))
    function _resolve_indices(ds::NCDataset, dnamesT::Tuple, roles::AbstractVector{Symbol}, Vraw)
        dnames = collect(String.(dnamesT))
        return _resolve_indices(ds, dnames, roles, Vraw)
    end
    function _resolve_indices(ds::NCDataset, dnames::AbstractVector{<:AbstractString}, roles::AbstractVector{Symbol}, Vraw)
        ix = findfirst(==( :lon  ), roles)
        iy = findfirst(==( :lat  ), roles)
        iz = findfirst(==( :z    ), roles)
        it = findfirst(==( :time ), roles)
        if ix === nothing
            ix = findfirst(n -> occursin(r"^x0*1$", lowercase(n)) || occursin(r"^x_?01$", lowercase(n)), dnames)
        end
        if iy === nothing
            iy = findfirst(n -> occursin(r"^x0*2$", lowercase(n)) || occursin(r"^x_?02$", lowercase(n)), dnames)
        end
        if iz === nothing
            iz = findfirst(n -> occursin(r"^x0*3$", lowercase(n)) || occursin(r"^x_?03$", lowercase(n)), dnames)
        end
        return ix, iy, iz, it
    end
    function coord_for_dim(ds::NCDataset, dimname::String, aliases::Vector{String}, size_in_dim::Int)
        if haskey(ds, dimname)
            return dimname, collect(ds[dimname][:])
        end
        for nm in aliases
            if haskey(ds, nm)
                dn = NCDatasets.dimnames(ds[nm])
                if !isempty(dn) && String(dn[1]) == dimname
                    return nm, collect(ds[nm][:])
                end
            end
        end
        return nothing, collect(1:size_in_dim)
    end
    if nd == 4
        ix, iy, iz, it = _resolve_indices(ds, dnames, roles, Vraw)
        ix === nothing && error("lon dim missing")
        iy === nothing && error("lat dim missing")
        iz === nothing && error("z dim missing")
        it === nothing && error("time dim missing")
        dim_lon = dnames[ix]; dim_lat = dnames[iy]; dim_z = dnames[iz]; dim_t = dnames[it]
        lonname, lon = coord_for_dim(ds, dim_lon, ["lon","longitude","x","grid_xt"], size(Vraw, ix))
        latname, lat = coord_for_dim(ds, dim_lat, ["lat","latitude","y","grid_yt"], size(Vraw, iy))
        zname  , z   = coord_for_dim(ds, dim_z,   ["alt","height","z","lev","level"], size(Vraw, iz))
        tname  = dim_t
        t = haskey(ds, dim_t) ? collect(ds[dim_t][:]) : collect(1:size(Vraw, it))
        perm = (ix, iy, iz, it)
        V    = perm == (1,2,3,4) ? Vraw : Array(PermutedDimsArray(Vraw, perm))
        return lat, lon, z, t, V, (latname === nothing ? dim_lat : latname,
                                   lonname === nothing ? dim_lon : lonname,
                                   zname   === nothing ? dim_z   : zname,
                                   tname)
    elseif nd == 3
        file_time === nothing && error("3D var needs file_time to synthesise time axis")
        ix, iy, iz, _  = _resolve_indices(ds, dnames, roles, Vraw)
        ix === nothing && error("lon dim missing")
        iy === nothing && error("lat dim missing")
        iz === nothing && error("z dim missing")
        dim_lon = dnames[ix]; dim_lat = dnames[iy]; dim_z = dnames[iz]
        lonname, lon = coord_for_dim(ds, dim_lon, ["lon","longitude","x","grid_xt"], size(Vraw, ix))
        latname, lat = coord_for_dim(ds, dim_lat, ["lat","latitude","y","grid_yt"], size(Vraw, iy))
        zname  , z   = coord_for_dim(ds, dim_z,   ["alt","height","z","lev","level"], size(Vraw, iz))
        V3 = (ix,iy,iz) == (1,2,3) ? Vraw : Array(PermutedDimsArray(Vraw, (ix,iy,iz)))
        V  = reshape(V3, size(V3,1), size(V3,2), size(V3,3), 1)
        t  = [file_time]
        return lat, lon, z, t, V, (latname === nothing ? dim_lat : latname,
                                   lonname === nothing ? dim_lon : lonname,
                                   zname   === nothing ? dim_z   : zname,
                                   "time")
    elseif nd == 2
        file_time === nothing && error("2D var needs file_time to synthesise time axis")
        ix, iy, _, _ = _resolve_indices(ds, dnames, roles, Vraw)
        ix === nothing && error("lon dim missing")
        iy === nothing && error("lat dim missing")
        dim_lon = dnames[ix]; dim_lat = dnames[iy]
        lonname, lon = coord_for_dim(ds, dim_lon, ["lon","longitude","x","grid_xt"], size(Vraw, ix))
        latname, lat = coord_for_dim(ds, dim_lat, ["lat","latitude","y","grid_yt"], size(Vraw, iy))
        V2 = (ix,iy) == (1,2) ? Vraw : Array(PermutedDimsArray(Vraw, (ix,iy)))
        V  = reshape(V2, size(V2,1), size(V2,2), 1, 1)
        zname = haskey(ds,"alt") ? "alt" : "level"
        z = [0.0]
        t = [file_time]
        return lat, lon, z, t, V, (latname === nothing ? dim_lat : latname,
                                   lonname === nothing ? dim_lon : lonname,
                                   zname, "time")
    else
        error("expected 3D or 4D variable, got ndims=$nd")
    end
end

# ============================ Time and altitude ==============================

function _decode_time_units(ds::NCDataset, tname::String, t::AbstractVector)
    if !haskey(ds, tname)
        return collect(t), nothing, nothing
    end
    units = get(ds[tname].attrib, "units", "")
    m = match(r"(seconds|minutes|hours|days) since (\d{4}-\d{2}-\d{2})(?:[ T](\d{2}:\d{2}:\d{2}))?", units)
    if m === nothing
        return collect(t), nothing, nothing
    end
    scale = m.captures[1]
    epoch_date = Date(m.captures[2])
    epoch_time = m.captures[3] === nothing ? Time(0) : Time(m.captures[3])
    epoch = DateTime(epoch_date, epoch_time)
    if eltype(t) <: DateTime
        tn = [_encode_query_time(tt, epoch, scale) for tt in t]
        return tn, epoch, scale
    else
        return collect(t), epoch, scale
    end
end

function _encode_query_time(dtq::DateTime, epoch::Union{DateTime,Nothing}, scale::Union{AbstractString,Nothing})
    epoch === nothing && return float(dtq.value)
    Δms = Dates.value(dtq - epoch)
    s = scale === nothing ? "days" : lowercase(String(scale))
    startswith(s,"sec") && return Δms/1_000
    startswith(s,"min") && return Δms/60_000
    startswith(s,"hour") && return Δms/3_600_000
    Δms/86_400_000
end

function _maybe_convert_alt(z::AbstractVector, alt_km::Real, ds::NCDataset, zname::String)
    units = get(ds[zname].attrib, "units", "")
    kind  = _classify_vertical_units(units)
    if kind === :km; return alt_km end
    if kind === :m;  return alt_km*1000 end
    if kind === :pressure; return :pressure end
    if kind === :index; error("vertical axis is index/level; pass a height variable to convert") end
    if kind === :missing || kind === :unknown; error("unsupported vertical units '$units'") end
end

# ====================== Interpolation (lon/lat/z/time) =======================

_nearest_index(vec::AbstractVector, x::Real) = findmin(abs.(vec .- x))[2]

function _bilinear_lonlat(lat::AbstractVector, lon::AbstractVector,
                          grid::AbstractArray{<:Real,2}, latq::Real, lonq::Real)
    lonq2 = _wrap_lon_for_grid(lon, lonq)
    ilat = clamp(searchsortedlast(lat, latq), 1, length(lat)-1)
    ilon = clamp(searchsortedlast(lon, lonq2), 1, length(lon)-1)
    φ1, φ2 = lat[ilat], lat[ilat+1]
    λ1, λ2 = lon[ilon], lon[ilon+1]
    tφ = (latq - φ1)/(φ2 - φ1); tλ = (lonq2 - λ1)/(λ2 - λ1)
    v11 = grid[ilon,   ilat  ]
    v21 = grid[ilon+1, ilat  ]
    v12 = grid[ilon,   ilat+1]
    v22 = grid[ilon+1, ilat+1]
    (1-tλ)*(1-tφ)*v11 + tλ*(1-tφ)*v21 + (1-tλ)*tφ*v12 + tλ*tφ*v22
end

function _interp3_linear(lat, lon, z, Vt, latq, lonq, zq)
    lonq2 = _wrap_lon_for_grid(lon, lonq)
    if length(z)==1
        Vz = Vt[:,:,1]
    else
        iz = clamp(searchsortedlast(z, zq), 1, length(z)-1)
        t = (zq - z[iz])/(z[iz+1]-z[iz])
        Vz = (1-t).*Vt[:,:,iz] .+ t.*Vt[:,:,iz+1]
    end
    _bilinear_lonlat(lat, lon, Vz, latq, lonq2)
end

function _interp3_logz_linear(lat, lon, z, Vt, latq, lonq, zq)
    lonq2 = _wrap_lon_for_grid(lon, lonq)
    if length(z) == 1
        return _bilinear_lonlat(lat, lon, Vt[:,:,1], latq, lonq2)
    end
    iz = clamp(searchsortedlast(z, zq), 1, length(z)-1)
    z1, z2 = z[iz], z[iz+1]
    if z1<=0 || z2<=0
        return _interp3_linear(lat,lon,z,Vt,latq,lonq2,zq)
    end
    V1 = Vt[:,:,iz]; V2 = Vt[:,:,iz+1]
    if any(x->x<=0 || !isfinite(x), V1) || any(x->x<=0 || !isfinite(x), V2)
        return _interp3_linear(lat,lon,z,Vt,latq,lonq2,zq)
    end
    t = (log(zq)-log(z1))/(log(z2)-log(z1))
    Vz = exp.((1-t).*log.(V1) .+ t.*log.(V2))
    _bilinear_lonlat(lat, lon, Vz, latq, lonq2)
end

function _quad_logz(z::AbstractVector, v::AbstractVector, zq::Real)
    mask = (z .> 0) .& isfinite.(z) .& (v .> 0) .& isfinite.(v)
    z2 = z[mask]; v2 = v[mask]
    if length(z2)==0; return NaN end
    if length(z2)==1; return v2[1] end
    x = log.(z2); y = log.(v2); xq = log(zq)
    ord = sortperm(abs.(x .- xq))
    sel = x[ord[1:min(3,end)]]; sely = y[ord[1:min(3,end)]]
    if length(sel)==2
        t = (xq-sel[1])/(sel[2]-sel[1]); return exp((1-t)*sely[1] + t*sely[2])
    end
    p = sortperm(sel); x1,x2,x3 = sel[p]; y1,y2,y3 = sely[p]
    d1=(x1-x2)*(x1-x3); d2=(x2-x1)*(x2-x3); d3=(x3-x1)*(x3-x2)
    if d1==0 || d2==0 || d3==0
        t = (xq - x1)/(x2 - x1); return exp((1-t)*y1 + t*y2)
    end
    L1=((xq-x2)*(xq-x3))/d1; L2=((xq-x1)*(xq-x3))/d2; L3=((xq-x1)*(xq-x2))/d3
    exp(L1*y1 + L2*y2 + L3*y3)
end

function _interp3_bilin_then_quadlogz(lat, lon, z, Vt, latq, lonq, zq)
    if length(z) <= 1
        return _bilinear_lonlat(lat, lon, Vt[:,:,1], latq, _wrap_lon_for_grid(lon, lonq))
    end
    vals = Vector{Float64}(undef, length(z))
    @inbounds for k in eachindex(z)
        vsl = view(Vt, :, :, k)
        vals[k] = _bilinear_lonlat(lat, lon, vsl, latq, lonq)
    end
    _quad_logz(z, vals, zq)
end

function _interp4(lat, lon, z, t, V, latq, lonq, zq, tq; mode::Symbol)
    lonq2 = _wrap_lon_for_grid(lon, lonq)
    if length(t)==1
        if mode == :linear
            return _interp3_linear(lat,lon,z,V[:,:,:,1],latq,lonq2,zq)
        elseif mode == :logz_linear
            return _interp3_logz_linear(lat,lon,z,V[:,:,:,1],latq,lonq2,zq)
        elseif mode == :logz_quadratic
            return _interp3_bilin_then_quadlogz(lat,lon,z,V[:,:,:,1],latq,lonq2,zq)
        else
            ilat=_nearest_index(lat,latq); ilon=_nearest_index(lon,lonq2); iz=_nearest_index(z,zq)
            return V[ilon,ilat,iz,1]
        end
    end
    if mode == :nearest
        ilat=_nearest_index(lat,latq); ilon=_nearest_index(lon,lonq2); iz=_nearest_index(z,zq); it=_nearest_index(t,tq)
        return V[ilon,ilat,iz,it]
    end
    it = clamp(searchsortedlast(t, tq), 1, length(t)-1)
    θ = (tq - t[it])/(t[it+1]-t[it])
    V1=V[:,:,:,it]; V2=V[:,:,:,it+1]
    v1 = mode==:linear ? _interp3_linear(lat,lon,z,V1,latq,lonq2,zq) :
         mode==:logz_linear ? _interp3_logz_linear(lat,lon,z,V1,latq,lonq2,zq) :
         _interp3_bilin_then_quadlogz(lat,lon,z,V1,latq,lonq2,zq)
    v2 = mode==:linear ? _interp3_linear(lat,lon,z,V2,latq,lonq2,zq) :
         mode==:logz_linear ? _interp3_logz_linear(lat,lon,z,V2,latq,lonq2,zq) :
         _interp3_bilin_then_quadlogz(lat,lon,z,V2,latq,lonq2,zq)
    (1-θ)*v1 + θ*v2
end

# ================= Pressure-level altitude mapping (gsm10 etc.) ==============

function _find_height_var(ds::NCDataset, zname::String)
    cands = String[]
    for k in keys(ds)
        v = ds[String(k)]
        nd = ndims(v)
        if nd in (3,4)
            nm = lowercase(String(k))
            a = Dict(v.attrib)
            units = lowercase(string(get(a,"units","")))
            std = lowercase(string(get(a,"standard_name","")))
            if (occursin("height", nm) || std=="height" || occursin("geopotential_height", std)) &&
               (occursin("m", units) || occursin("meter", units))
                push!(cands, String(k))
            end
        end
    end
    isempty(cands) && error("no height variable found to map altitude on pressure levels")
    first(cands)
end

function _alt_to_zindex(ds::NCDataset, lat, lon, z, t, V, latq, lonq, alt_m, zname::String)
    hvar = _find_height_var(ds, zname)
    H = ds[hvar]
    Hraw = Array(H)
    dnames = String.(NCDatasets.dimnames(H))
    roles = map(n -> lowercase(n)=="time" ? :time :
                    occursin("lon",lowercase(n)) || n in ("x","grid_xt","i") ? :lon :
                    occursin("lat",lowercase(n)) || n in ("y","grid_yt","j") ? :lat :
                    :z, dnames)
    perm = (findfirst(==( :lon ), roles),
            findfirst(==( :lat ), roles),
            findfirst(==( :z   ), roles),
            findfirst(==( :time), roles))
    if any(x->x===nothing, perm); error("height variable dims incompatible") end
    Hn = ndims(Hraw)==4 ? Array(PermutedDimsArray(Hraw, Tuple(perm))) :
         ndims(Hraw)==3 ? reshape(Array(PermutedDimsArray(Hraw, Tuple(perm[1:3]))), size(Hraw,perm[1]), size(Hraw,perm[2]), size(Hraw,perm[3]), 1) :
         error("height variable ndims must be 3 or 4")
    it = 1
    ilat=_nearest_index(lat,latq); ilon=_nearest_index(lon,_wrap_lon_for_grid(lon,lonq))
    h_prof = vec(view(Hn, ilon, ilat, :, it))
    iz = clamp(searchsortedlast(h_prof, alt_m), 1, length(h_prof)-1)
    tζ = (alt_m - h_prof[iz])/(h_prof[iz+1]-h_prof[iz])
    (iz, tζ)
end

# ================================ Validation ================================

function _validate_query_args(itp::WFSInterpolator, dt::DateTime, latq::Real, lonq::Real, alt_km::Real)
    mode = _normalize_interp(itp.interpolation)
    mode in _ALLOWED_INTERP_NORM || error("interpolation must be one of $(collect(_ALLOWED_INTERP_NORM)) or :sciml")
    isfinite(latq) && -90 <= latq <= 90 || error("lat out of range")
    isfinite(lonq) || error("lon not finite")
    isfinite(alt_km) || error("alt_km not finite")
    (alt_km > 0) || error("alt_km must be > 0")
    mode
end

# ================================== API ======================================

"""
    get_value(itp::WFSInterpolator, dt::DateTime, lon::Real, lat::Real, alt_km::Real; varname=itp.varname) -> Float64

Interpolated value of `varname` from the chosen product/stream at `(dt, lon, lat, alt_km)`.
If the vertical axis is on pressure levels, a height field is used to map altitude to level.
Temporal interpolation is linear between the two nearest valid files.
"""
function get_value(itp::WFSInterpolator, dt::DateTime, lonq::Real, latq::Real, alt_km::Real; varname::String=itp.varname)
    mode = _validate_query_args(itp, dt, latq, lonq, alt_km)
    url_lo, url_hi, t_lo, t_hi = _pick_two_files(itp, dt)
    path_lo = _download_http_cached(itp, url_lo)
    path_hi = _download_http_cached(itp, url_hi)
    ds_lo = NCDataset(path_lo, "r"); ds_hi = NCDataset(path_hi, "r")
    try
        lat, lon, z, t, V, (latname, lonname, zname, tname) = _load_grids(ds_lo, varname; file_time=t_lo)
        tdts, epoch, scale = _decode_time_units(ds_lo, tname, t)
        tq_lo = epoch === nothing ? t_lo : _encode_query_time(t_lo, epoch, scale)
        alt_q = _maybe_convert_alt(z, alt_km, ds_lo, zname)
        v_lo = if alt_q === :pressure
            iz, τ = _alt_to_zindex(ds_lo, lat, lon, z, tdts, V, latq, lonq, alt_km*1000, zname)
            v_iz  = _interp4(lat, lon, z, tdts, V, latq, lonq, z[iz],   tq_lo; mode=:linear)
            v_iz1 = _interp4(lat, lon, z, tdts, V, latq, lonq, z[iz+1], tq_lo; mode=:linear)
            (1-τ)*v_iz + τ*v_iz1
        else
            _interp4(lat, lon, z, tdts, V, latq, lonq, alt_q, tq_lo; mode=mode)
        end
        lat2, lon2, z2, t2, V2, (latname2, lonname2, zname2, tname2) = _load_grids(ds_hi, varname; file_time=t_hi)
        tdts2, epoch2, scale2 = _decode_time_units(ds_hi, tname2, t2)
        tq_hi = epoch2 === nothing ? t_hi : _encode_query_time(t_hi, epoch2, scale2)
        alt_q2 = _maybe_convert_alt(z2, alt_km, ds_hi, zname2)
        v_hi = if alt_q2 === :pressure
            iz2, τ2 = _alt_to_zindex(ds_hi, lat2, lon2, z2, tdts2, V2, latq, lonq, alt_km*1000, zname2)
            v_iz  = _interp4(lat2, lon2, z2, tdts2, V2, latq, lonq, z2[iz2],   tq_hi; mode=:linear)
            v_iz1 = _interp4(lat2, lon2, z2, tdts2, V2, latq, lonq, z2[iz2+1], tq_hi; mode=:linear)
            (1-τ2)*v_iz + τ2*v_iz1
        else
            _interp4(lat2, lon2, z2, tdts2, V2, latq, lonq, alt_q2, tq_hi; mode=mode)
        end
        if t_lo == t_hi; return float(v_lo) end
        itp_t = DataInterpolations.LinearInterpolation([float(v_lo), float(v_hi)],
                                                       [Dates.value(t_lo), Dates.value(t_hi)])
        return itp_t(Dates.value(dt))
    finally
        close(ds_lo); close(ds_hi)
    end
end

"""
    get_batch(itp, dts, lons, lats, alts_km; varname=itp.varname) -> Vector{Float64}

Vectorised wrapper over [`get_value`](@ref). All inputs must have the same length.
"""
function get_batch(itp::WFSInterpolator, dts::AbstractVector{<:DateTime},
                   lons::AbstractVector, lats::AbstractVector, alts_km::AbstractVector; varname::String=itp.varname)
    n = length(dts)
    @assert length(lons)==n==length(lats)==length(alts_km)
    [get_value(itp, dts[i], lons[i], lats[i], alts_km[i]; varname=varname) for i in 1:n]
end

"""
    list_vars(itp::WFSInterpolator, dt::DateTime) -> Vector{String}

List 3-D variables available in the low-bracketing file for the given `dt`.
"""
function list_vars(itp::WFSInterpolator, dt::DateTime)
    url_lo, _, _, _ = _pick_two_files(itp, dt)
    p_lo = _download_http_cached(itp, url_lo)
    ds = NCDatasets.NCDataset(p_lo, "r")
    try
        v3 = String[]
        for k in keys(ds)
            v = ds[String(k)]
            ndims(v) == 3 && push!(v3, String(k))
        end
        sort(v3)
    finally
        close(ds)
    end
end

"""
    dump_sample(itp, dt, lonq, latq, alt_km;
                modes = (:nearest, :linear, :logz_linear, :logz_quadratic, :sciml),
                varfilter = :all) -> Dict

Enumerate 3-D variables at the query point and report per-file values for each interpolation mode,
plus the linear time-blended value. `varfilter` may be `:all`, a `Vector{String}`, or a `Regex`.
Returns a dictionary containing file paths, times, grid metadata, and values.
"""
function dump_sample(itp::WFSInterpolator, dt::DateTime, lonq::Real, latq::Real, alt_km::Real;
                     modes = (:nearest, :linear, :logz_linear, :logz_quadratic, :sciml),
                     varfilter = :all)
    url_lo, url_hi, t_lo, t_hi = _pick_two_files(itp, dt)
    path_lo = _download_http_cached(itp, url_lo)
    path_hi = _download_http_cached(itp, url_hi)
    ds_lo = NCDataset(path_lo, "r"); ds_hi = NCDataset(path_hi, "r")
    try
        function three_d_names(ds::NCDataset)
            out = String[]
            for k in keys(ds)
                v = ds[String(k)]
                nd = ndims(v)
                if nd == 2 || nd == 3
                    push!(out, String(k))
                end
            end
            out
        end
        cand = intersect(Set(three_d_names(ds_lo)), Set(three_d_names(ds_hi))) |> collect |> sort
        if varfilter isa Vector{String}
            cand = [v for v in cand if v in varfilter]
        elseif varfilter isa Regex
            cand = [v for v in cand if occursin(varfilter, v)]
        elseif varfilter !== :all
            error("varfilter must be :all, Vector{String}, or Regex")
        end
        isempty(cand) && error("no 3-D variables found to dump")
        latL, lonL, zL, tL, VL, (latname, lonname, zname, tname) =
            _load_grids(ds_lo, cand[1]; file_time=t_lo)
        tdtsL, epochL, scaleL = _decode_time_units(ds_lo, tname, tL)
        tq_lo = epochL === nothing ? t_lo : _encode_query_time(t_lo, epochL, scaleL)
        alt_axis_lo = _maybe_convert_alt(zL, alt_km, ds_lo, zname)
        lonq2 = _wrap_lon_for_grid(lonL, lonq)
        ilon = _nearest_index(lonL, lonq2)
        ilat = _nearest_index(latL, latq)
        iz   = alt_axis_lo === :pressure ? missing : _nearest_index(zL, alt_axis_lo)
        latH, lonH, zH, tH, VH, (latnameH, lonnameH, znameH, tnameH) =
            _load_grids(ds_hi, cand[1]; file_time=t_hi)
        tdtsH, epochH, scaleH = _decode_time_units(ds_hi, tnameH, tH)
        tq_hi = epochH === nothing ? t_hi : _encode_query_time(t_hi, epochH, scaleH)
        alt_axis_hi = _maybe_convert_alt(zH, alt_km, ds_hi, znameH)
        w = (Dates.value(dt) - Dates.value(t_lo)) / (Dates.value(t_hi) - Dates.value(t_lo))
        results = Dict{String,Any}()
        for var in cand
            res_var = Dict{String,Any}()
            latL, lonL, zL, tL, VL, (latname, lonname, zname, tname) =
                _load_grids(ds_lo, var; file_time=t_lo)
            tdtsL, epochL, scaleL = _decode_time_units(ds_lo, tname, tL)
            tq_lo = epochL === nothing ? t_lo : _encode_query_time(t_lo, epochL, scaleL)
            alt_q_lo = _maybe_convert_alt(zL, alt_km, ds_lo, zname)
            latH, lonH, zH, tH, VH, (latnameH, lonnameH, znameH, tnameH) =
                _load_grids(ds_hi, var; file_time=t_hi)
            tdtsH, epochH, scaleH = _decode_time_units(ds_hi, tnameH, tH)
            tq_hi = epochH === nothing ? t_hi : _encode_query_time(t_hi, epochH, scaleH)
            alt_q_hi = _maybe_convert_alt(zH, alt_km, ds_hi, znameH)
            for mode in modes
                modeN = _normalize_interp(mode)
                v_lo = if alt_q_lo === :pressure
                    iz0, τ0 = _alt_to_zindex(ds_lo, latL, lonL, zL, tdtsL, VL, latq, lonq, alt_km*1000, zname)
                    v_iz  = _interp4(latL, lonL, zL, tdtsL, VL, latq, lonq, zL[iz0],   tq_lo; mode=:linear)
                    v_iz1 = _interp4(latL, lonL, zL, tdtsL, VL, latq, lonq, zL[iz0+1], tq_lo; mode=:linear)
                    (1-τ0)*v_iz + τ0*v_iz1
                else
                    _interp4(latL, lonL, zL, tdtsL, VL, latq, lonq, alt_q_lo, tq_lo; mode=modeN)
                end
                v_hi = if alt_q_hi === :pressure
                    iz1, τ1 = _alt_to_zindex(ds_hi, latH, lonH, zH, tdtsH, VH, latq, lonq, alt_km*1000, znameH)
                    v_iz  = _interp4(latH, lonH, zH, tdtsH, VH, latq, lonq, zH[iz1],   tq_hi; mode=:linear)
                    v_iz1 = _interp4(latH, lonH, zH, tdtsH, VH, latq, lonq, zH[iz1+1], tq_hi; mode=:linear)
                    (1-τ1)*v_iz + τ1*v_iz1
                else
                    _interp4(latH, lonH, zH, tdtsH, VH, latq, lonq, alt_q_hi, tq_hi; mode=modeN)
                end
                v_blend = (1-w)*float(v_lo) + w*float(v_hi)
                res_var[string(modeN)] = Dict("lo"=>float(v_lo), "hi"=>float(v_hi), "blend"=>float(v_blend))
            end
            units = get(ds_lo[var].attrib, "units", "")
            res_var["units"] = String(units)
            results[var] = res_var
        end
        return Dict(
            "files" => Dict("lo"=>path_lo, "hi"=>path_hi),
            "times" => Dict("lo"=>t_lo, "hi"=>t_hi, "w"=>w),
            "grid"  => Dict("lon"=>lonL, "lat"=>latL, "z"=>zL,
                            "nearest"=>Dict("ilon"=>ilon, "ilat"=>ilat, "iz"=>iz, "lonq_mapped"=>lonq2)),
            "vars"  => results
        )
    finally
        close(ds_lo); close(ds_hi)
    end
end

"""
    dump_all(dt, lon, lat, alt_km;
             products = ["wfs","wrs"],
             streams  = ["ipe05","ipe10","gsm05","gsm10"],
             modes    = (:nearest, :linear, :logz_linear, :logz_quadratic, :sciml),
             varfilter = :all,
             cache_dir = "cache") -> Dict

Run [`dump_sample`](@ref) across multiple product/stream combinations and collect results.
Errors per combination are captured and returned in the output dictionary.
"""
function dump_all(dt::DateTime, lon::Real, lat::Real, alt_km::Real;
                  products = ["wfs","wrs"],
                  streams  = ["ipe05","ipe10","gsm05","gsm10"],
                  modes    = (:nearest, :linear, :logz_linear, :logz_quadratic, :sciml),
                  varfilter = :all,
                  cache_dir = "cache")
    results = Dict{String,Any}()
    for prod in products
        for strm in streams
            key = "$(prod):$(strm)"
            itp = WFSInterpolator(product=prod, stream=strm, cache_dir=cache_dir)
            try
                rep = dump_sample(itp, dt, lon, lat, alt_km; modes=modes, varfilter=varfilter)
                results[key] = rep
            catch err
                @warn "Failed for $key" error=err
                results[key] = Dict("error" => string(err))
            end
        end
    end
    return results
end

end # module

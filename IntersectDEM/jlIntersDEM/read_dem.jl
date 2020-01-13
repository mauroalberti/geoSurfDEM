

"""
paths
"""

input_project_fldrpth = "/home/mauro-grigio/Documents/projects/geoSurfDEM" 

"""
functions
"""

mutable struct AsciiDemParams
    ncols::Int32
    nrows::Int32
    xllcorner::Float64
    yllcorner::Float64
    cellsize::Float64
    nodata_value::Float64
end


function parse_line(file)

	return split(readline(file))[2]

end


function read_header(fltpth)

	ncols::Int32 = 0
	nrows::Int32 = 0
	xllcorner::Float64 = 0.0
	yllcorner::Float64 = 0.0
	cellsize::Float64 = 0.0
	nodata_value::Float64 = 0.0

	open(dem_flpth) do fdem

		ncols = tryparse(Int32, parse_line(fdem))
		nrows = tryparse(Int32, parse_line(fdem))
		xllcorner = tryparse(Float64, parse_line(fdem))
		yllcorner = tryparse(Float64, parse_line(fdem))
		cellsize = tryparse(Float64, parse_line(fdem))
		nodata_value = tryparse(Float64, parse_line(fdem))

	end 

	return AsciiDemParams(
		ncols,
		nrows,
		xllcorner,
		yllcorner,
		cellsize,
		nodata_value
	)

end


function progr2array(cntr, ncols)

	curr_row, curr_col = divrem(cntr - 1, ncols)
	curr_row += 1
	curr_col += 1

	return curr_row, curr_col


function read_dem_values(fltpth, dem_params::AsciiDemParams)

	dem_matrix = Matrix{Float64}(undef, dem_params.nrows, dem_params.ncols)

	open(dem_flpth) do fdem

	    linecounter = 0
	    nvals = 0
	    timetaken = @elapsed for l in eachline(fdem)
	        linecounter += 1
	        #println("\t", linecounter, "\n", l)
	        if linecounter < 7
	        	continue
	        end
	        vals = split(l)
	        for val in vals
	        	nvals += 1
	        	curr_row, curr_col = progr2array(nvals, dem_params.ncols))
	        	println(nvacurr_row, curr_col = ls, curr_row, curr_col)
	        end
	    end

	end

	#println(totaltime, "\n\t", totallines)

	#return dem_matrix

end


data_fldrpth = "test_data/IntersectDEM"
dem_flnm = "dem_malpi_aster_wgs84utm33.asc"

dem_flpth = input_project_fldrpth * "/" * data_fldrpth * "/" * dem_flnm

# read dem header

dem_header = read_header(dem_flpth)

#print(dem_header)

read_dem_values(dem_flpth, dem_header)

#print(dem_data)


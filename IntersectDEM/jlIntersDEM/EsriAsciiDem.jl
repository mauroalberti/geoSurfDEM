
module EsriAsciiDem

	export read_header, read_dem_values

	const num_header_rows = 7


	struct AsciiDemParams
	    ncols::Int32
	    nrows::Int32
	    xllcorner::Float64
	    yllcorner::Float64
	    cellsize::Float64
	    nodata_value::String
	end


	function parse_line(file::String) :: String

		return split(readline(file))[2]

	end


	function read_header(fltpth::String) :: AsciiDemParams

		ncols::Int32 = 0
		nrows::Int32 = 0
		xllcorner::Float64 = 0.0
		yllcorner::Float64 = 0.0
		cellsize::Float64 = 0.0
		nodata_value::String = ""

		open(fltpth) do fdem

			ncols = tryparse(Int32, parse_line(fdem))
			nrows = tryparse(Int32, parse_line(fdem))
			xllcorner = tryparse(Float64, parse_line(fdem))
			yllcorner = tryparse(Float64, parse_line(fdem))
			cellsize = tryparse(Float64, parse_line(fdem))
			nodata_value = parse_line(fdem)

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


	function read_values(fltpth::String, dem_params::AsciiDemParams) :: Matrix{Float64}

		dem_matrix = Matrix{Float64}(undef, dem_params.nrows, dem_params.ncols)

		open(fltpth) do fdem

		    linecounter = 0
		    nvals = 0
		    for l in eachline(fdem)
		        linecounter += 1
		        if linecounter < num_header_rows
		        	continue
		        end
		        vals = split(l)
		        for val in vals
		        	nvals += 1
		        	curr_row, curr_col = progr2array(nvals, dem_params.ncols))		        	
		        	parsed_val = val == dem_params.nodata_value ? missing : Float64(val) 
		        	dem_matrix[curr_row, curr_col] = parsed_val
		        end
		    end

		end

		return dem_matrix

	end


end


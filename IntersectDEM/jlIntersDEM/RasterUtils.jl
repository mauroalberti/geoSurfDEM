
module RasterUtils


	function progr2array(cntr::Int32, ncols::Int32) :: Tuple{Int32,Int32}

		curr_row, curr_col = divrem(cntr - 1, ncols)
		curr_row += 1
		curr_col += 1

		return curr_row, curr_col


end
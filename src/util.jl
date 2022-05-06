using DataFrames
import Glob 

## Use abstract type to select range conditions
## Method inspired by https://www.juliabloggers.com/julia-dispatching-enum-versus-type/
abstract type ValueBound end
struct OpenBound   <: ValueBound end
struct ClosedBound <: ValueBound end
struct LeftClosed  <: ValueBound end
struct RightClosed <: ValueBound end


function getbolddirs(bdir::AbstractString)
	fdrs = Glob.glob("*", bdir)
	[split(f,"/")[end] for f in fdrs]
end


function findpattern(nxfiles::Vector{String}, pattern::String, spl="_", pos=1)
	REPS = []
	for f in nxfiles
		fsx = split(f, spl)
		fa = findall(x->x==pattern, fsx)
		if length(fa) > 0
			indx = fa[1] + pos
			push!(REPS, fsx[indx])
		end
	end
	unique(REPS)
end

function range_bound(xmin::Real, xmax::Real, ::Type{OpenBound  })
	x -> (x >  xmin) & (x <  xmax)
end

function range_bound(xmin::Real, xmax::Real, ::Type{ClosedBound})
	x -> (x >= xmin) & (x <= xmax)
end

function range_bound(xmin::Real, xmax::Real, ::Type{LeftClosed })
	x -> (x >= xmin) & (x <  xmax)
end

function range_bound(xmin::Real, xmax::Real, ::Type{RightClosed})
	x -> (x >  xmin) .& (x <= xmax)
end


"""
	swap(x1::T, x2::T, cond::Bool) where T
	swaps x and y if cond is false
"""
function swap(x1::T, x2::T, cond::Bool) where T
    if cond
        return x1, x2
    else
        return x2, x1
    end
end

# Vector and data frames

"""
	in_range(x, xmin, xmax)

Given vector x, select values between xmin and xmax
"""
function in_range(x::Vector{<:Real}, xmin::Real, xmax::Real,
				  interval::Type{<:ValueBound}=OpenBound)
	mask = broadcast(range_bound(xmin, xmax, interval), x)
    return x[mask]
end


"""
	select_values

Generic function to get values in a DataFrame given a condition
function.

** Arguments **
dbdf     ::DataFrame
cond_func::Function Should take a value and return Boolean (brodcasting in function)
column   ::String The column to be used with cond_func
TODO: Change access from String to Symbol
TODO: Is it possible to specify that cond must take Vector and return Bool?
"""
function select_values(dbdf::DataFrame, cond_func::Function, column::String)
	mask = broadcast(cond_func, dbdf[!, column])
	return dbdf[mask, :]
end


"""
	select_event(dbdf::DataFrame, index::Int64)

Take the event dataframe and the index of an event and returns a data frame
which selects that particular event

"""
function select_event(dbdf::DataFrame, event_id::Integer)
	return select_by_column_value(dbdf, "event_id", event_id)
end


"""
	select_by_column_value(df::DataFrame, column::String, value)

Select elements in the DF which have "value" in "column"
"""
function select_by_column_value(df::DataFrame, column::String, value)
	return select_values(df, x -> x == value, column)
end

"""
	select_by_column_value_lt(df::DataFrame, column::String, value)

Select elements in the DF which are less than "value" in "column"
"""
function select_by_column_value_lt(df::DataFrame, column::String, value)
	return select_values(df, x -> x < value, column)
end


"""
	select_by_column_value_le(df::DataFrame, column::String, value)

Select elements in the DF which are less or equal than "value" in "column"
"""
function select_by_column_value_le(df::DataFrame, column::String, value)
	return select_values(df, x -> x <= value, column)
end


"""
	select_by_column_value_gt(df::DataFrame, column::String, value)

Select elements in the DF which are larger than "value" in "column"
"""
function select_by_column_value_gt(df::DataFrame, column::String, value)
	return select_values(df, x -> x > value, column)
end


"""
	select_by_column_value_ge(df::DataFrame, column::String, value)

Select elements in the DF which are larger or equal than "value" in "column"
"""
function select_by_column_value_ge(df::DataFrame, column::String, value)
	return select_values(df, x -> x >= value, column)
end


"""
	select_by_column_value_interval(df::DataFrame, column::String, valuef, valuel)

Select elements in the DF which are in interval (valuef, valuel)
"""
function select_by_column_value_interval(df::DataFrame, column::String, valuef, valuel)
	cond_func = range_bound(valuef, valuel, OpenBound)
	return select_values(df, cond_func, column)
end


"""
	select_by_column_value_closed_interval(df::DataFrame, column::String, valuef, valuel)

Select elements in the DF which are in interval [valuef, valuel]
"""
function select_by_column_value_closed_interval(df::DataFrame, column::String, valuef, valuel)
	cond_func = range_bound(valuef, valuel, ClosedBound)
	return select_values(df, cond_func, column)
end


"""
	select_by_column_value_closed_left_interval(df::DataFrame, column::String, valuef, valuel)

Select elements in the DF which are in interval [valuef, valuel)
"""
function select_by_column_value_closed_left_interval(df::DataFrame, column::String, valuef, valuel)
	cond_func = range_bound(valuef, valuel, LeftClosed)
	return select_values(df, cond_func, column)
end


"""
	select_by_column_value_closed_right_interval(df::DataFrame, column::String, valuef, valuel)

Select elements in the DF which are in interval (valuef, valuel]
"""
function select_by_column_value_closed_right_interval(df::DataFrame, column::String, valuef, valuel)
	cond_func = range_bound(valuef, valuel, RightClosed)
	return select_values(df, cond_func, column)
end



"""
	select_by_index(df::DataFrame, column::String, value::Integer)

Select elements in the DF which have "value" (Integer) in "column"
!Name is misleading, does not select by index.
"""
function select_by_index(df::DataFrame, column::String, value::Integer)
	return select_by_column_value(df, column, value)
end



""""
	find_max_xy(df, xc, yc)

Return ymax and x such that ymax = f(x).

**Description:**

In a DataFrame one has often "XY" variables, that is, a pair of columns
"X" and "Y" which represent correlated variables (e.g intensity and wavelength).
In such cases one often wants the XY maximum, that is, finding the maximum
of "Y" (call it ymax) and the corresponding value in "X"
(x for which y is ymax). This corresponds, for example, to the wavelength at
which the intensity is maximal.


**Arguments:**
- `df::DataFrame`: data frame holding the data.
- `xc::String`: the X column.
- `yc::String`: the Y column.
"""
find_max_xy(df::DataFrame, xc::String, yc::String) = find_max_xy(df[!, xc], df[!, yc])

function find_max_xy(x::Vector{<:Real}, y::Vector{<:Real})
	ymax, imax = findmax(y)
	x_ymax     = x[imax]
	return imax, x_ymax, ymax
end

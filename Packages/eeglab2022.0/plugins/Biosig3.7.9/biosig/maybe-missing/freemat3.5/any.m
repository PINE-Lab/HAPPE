%!
%@Module ANY Any True Function
%@@Section ELEMENTARY
%@@Usage
%Reduces a logical array along a given dimension by testing for any
%logical 1s.  The general 
%syntax for its use is
%@[
%  y = any(x,d)
%@]
%where @|x| is an @|n|-dimensions array of @|logical| type.
%The output is of @|logical| type.  The argument @|d| is 
%optional, and denotes the dimension along which to operate.
%The output @|y| is the same size as @|x|, except that it is 
%singular along the operated direction.  So, for example,
%if @|x| is a @|3 x 3 x 4| array, and we @|any| operation along
%dimension @|d=2|, then the output is of size @|3 x 1 x 4|.
%@@Function Internals
%The output is computed via
%\[
%y(m_1,\ldots,m_{d-1},1,m_{d+1},\ldots,m_{p}) = 
%\max_{k} x(m_1,\ldots,m_{d-1},k,m_{d+1},\ldots,m_{p})
%\]
%If @|d| is omitted, then the summation is taken along the 
%first non-singleton dimension of @|x|. 
%@@Example
%The following piece of code demonstrates various uses of the summation
%function
%@<
%A = [1,0,0;1,0,0;0,0,1]
%@>
%We start by calling @|any| without a dimension argument, in which 
%case it defaults to the first nonsingular dimension (in this case, 
%along the columns or @|d = 1|).
%@<
%any(A)
%@>
%Next, we apply the @|any| operation along the rows.
%@<
%any(A,2)
%@>
%@@Tests
%@$"y=any([1,0,0;1,0,0;0,0,1],2)","[1;1;1]","exact"
%!

% Copyright (c) 2002-2007 Samit Basu
% Copyright (c) 2008 Alois Schloegl
% Licensed under GPL v2 or later

function y = any(A,DIM)
	if (nargin == 0)
		error('any function requires at least one argument');
	elseif (nargin == 1)
		DIM = [];    
	end

	sz = size(A); 	
	if isempty(DIM),
		DIM=min(find(sz ~= 1));
		if (DIM<1), DIM = 1; end;
		if isempty(DIM), DIM = 1; end;
	end;

	if DIM>length(sz),
		sz = [sz,ones(1,DIM-length(sz))];
	end;

	%% handling of N-dim arrays
	D1 = prod(sz(1:DIM-1));
	D3 = prod(sz(DIM+1:length(sz)));
	y  = repmat(logical(0),[sz(1:DIM-1),1,sz(DIM+1:length(sz))]);
	A(isnan(A)) = 0;
	for k = 0:D1-1,
	for l = 0:D3-1,
	        xi = k + l * D1*sz(DIM) + 1 ;
		xo = k + l * D1 + 1;
	        f = logical(0); 
	        m = 0;
	        while (~f && m<sz(DIM))		% short-cut evaluation 
			f = f || (A(xi+D1*m)~=0);
			m = m+1;
		end;
		y(xo) = f;
	end;
	end;

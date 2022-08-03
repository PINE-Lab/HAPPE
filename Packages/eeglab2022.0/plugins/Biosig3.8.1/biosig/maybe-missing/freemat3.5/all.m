%!
%@Module ALL All True Function
%@@Section ELEMENTARY
%@@Usage
%Reduces a logical array along a given DIMension by testing for all
%logical 1s.  The general 
%syntax for its use is
%@[
%  y = all(x,d)
%@]
%where @|x| is an @|n|-DIMensions array of @|logical| type.
%The output is of @|logical| type.  The argument @|d| is 
%optional, and denotes the DIMension along which to operate.
%The output @|y| is the same size as @|x|, except that it is 
%singular along the operated direction.  So, for example,
%if @|x| is a @|3 x 3 x 4| array, and we @|all| operation along
%DIMension @|d=2|, then the output is of size @|3 x 1 x 4|.
%@@Function Internals
%The output is computed via
%\[
%y(m_1,\ldots,m_{d-1},1,m_{d+1},\ldots,m_{p}) = 
%\min_{k} x(m_1,\ldots,m_{d-1},k,m_{d+1},\ldots,m_{p})
%\]
%If @|d| is omitted, then the minimum is taken over all elements of
%@|x|.
%@@Example
%The following piece of code demonstrates various uses of the @|all|
%function
%@<
%A = [1,0,0;1,0,0;0,0,1]
%@>
%We start by calling @|all| without a DIMension argument, in which 
%case it defaults to testing all values of @|A|
%@<
%all(A)
%@>
%The @|all| function is useful in expressions also.
%@<
%all(A>=0)
%@>
%Next, we apply the @|all| operation along the rows.
%@<
%all(A,2)
%@>
%@@Tests
%@$"y=all([1,0,0;1,1,1;0,1,1],2)","[0;1;0]","exact"
%!

% Copyright (c) 2002-2007 Samit Basu
% Copyright (c) 2008 Alois Schloegl
% Licensed under GPL v2 or later

function y = all(A,DIM)
	if (nargin == 0)
		error('all function requires at least one argument');
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

	D1 = prod(sz(1:DIM-1));
	D3 = prod(sz(DIM+1:length(sz)));
	y  = repmat(logical(1),[sz(1:DIM-1),1,sz(DIM+1:length(sz))]);
	for k = 0:D1-1,
	for l = 0:D3-1,
	        xi = k + l * D1*sz(DIM) + 1 ;
		xo = k + l * D1 + 1;
	        f = logical(1); 
	        m = 0;
	        while (f && m<sz(DIM))
			f = f && (A(xi+D1*m)~=0);
			m = m+1;
		end;
		y(xo) = f;
	end;
	end;
	
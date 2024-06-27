function [TOPT] = optimum_isi_spike_burst_separation(ISI, ORDER);
%  OPTIMUM_ISI_SPIKE_BURST_SEPARATION tries to find an optimum
%  threshold to separate spike from burst intervals. The problem is
%  described in [1](pp.534-6). There are different proposed [2-5]
%  described by equ. (45)-(47) in [1].
%  All these methods are based on fitting a double exponential
%  function to a binned histogram. The method used here does not
%  require binned histograms; and is much more robust, especially
%  when small sample sizes are used.
%
%  The method is based on a polynomial fit of the log-logit transformed
%  cumulative distribution function (CDF) of a given set of inter-spike
%  intervals. The minimum of the first derivative is returned as
%  the optimum threshold. The 2nd derivative is used to exclude the
%  influnce of extremes/outliers.
%
%  Usage:
%    Topt = optimum_isi_spike_burst_separation(ISI);
%    Topt = optimum_isi_spike_burst_separation(ISI, ORDER);
%
%  ISI   list (i.e. vector) of interspike intervals
%  ORDER order of polynomial used to approximate the
%        linearized CDF function [default: 7]
%  Topt  optimum time interval for identifying burst intervals.
%
% see also:  DETECT_SPIKES_BURSTS, SPIKE2BURSTS
%
% References:
% [1] David Colquhoun and F.J. Sigworth (1995)
%   Chapter 19. Fitting and Statistical Analysis of Single-Channel Records.
%   Eds. Bert Sakmann and Erwin Neher "Single-Channel Recording", 2nd Ed.
%   Plenum Press, New York and London. p.535.
% [2] Jackson et al. (1983)
% [3] Magleby and Pallotta (1983)
% [4] Clapham and Neher (1984)
% [5] Colquhoun and Sakmann (1985)


%    Copyright (C) 2014 by Alois Schloegl <alois.schloegl@ist.ac.at>
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.



VERBOSE=0; % 0: no display, 1: display


if (nargin<1) || (ischar(ISI) && strcmpi(ISI,'demo'));

	%% generate some demo data - used for testing the method
	N1=1000;
	N2=100;
	M1=0;
	M2=3;
	S1=1;
	S2=1;
	x=exp([S1*randn(N1,1)+M1; S2*randn(N2,1)+M2]);
	if (VERBOSE > 0) 	
		plot(x,'x')
	end
	% The outcome for log(Topt) is in the range of 2.6 +- 0.22
else
	assert(all(ISI>0));
	x = ISI;
end;

if (nargin<2) || isempty(ORDER),
	ORDER = 7; %% order of the polynomial
end;


HIS=histo3(x);

if (VERBOSE > 0) 	
	subplot(121);
	HIS.N = sum(HIS.H);
	semilogx([1e-4;HIS.X;inf]', logit([0;cumsum(HIS.H);HIS.N]/HIS.N), 'b');			
	xlabel('inter-spike interval [s]');
	ylabel('logit(CDF) {blue}');
	title('');
	grid on;
end;

% transform X and Y to a linear scale		
Y = logit(cumsum(HIS.H)/HIS.N);
%Y = probit(cumsum(HIS.H)/HIS.N);


X = log(HIS.X);
T = [min(X):(max(X)-min(X))/max(1000, HIS.N*2):max(X)]';


% perform polynomial fit of log-logit transformed CDF
[P,S,M] = polyfit(X(1:end-1),Y(1:end-1), ORDER);


% 1st, 2nd and 3rd derivative polynomial
P1 = polyder(P);
P2 = polyder(P1);
% P3 = polyder(P2);


% approximate 1st derivative function
a = polyval(P1,T,[],M);

% while the 2nd derivative must be close to zero - threshold of 1.0 seems to be good enough for excluding extreme outcomes
a ( abs(polyval(P2,T,[],M) ) > 1.0 ) = NaN;

% while the 3rd derivative must be smaller than 0, otherwise its not a minimum
% a ( abs(polyval(P3,T,[],M) ) <= 0 ) = NaN;


% find minimum of the 1st derivative
% 2nd derivative must contain a zero crossing (must be close to zero)
[tmp,ix] = min(a);

% transform to original scale on x-axis. 	
TOPT = exp(T(ix));

if (VERBOSE > 0)
	hold on
	plot(TOPT*[1;1],[-5;5],'x-k')
	hold off

	subplot(122);
		
	semilogx(exp(X),[Y,polyval(P,X,[],M),polyval(P1,X,[],M),polyval(P2,X,[],M)],'-',exp(T),[polyval(P,T,[],M),polyval(P1,T,[],M),polyval(P2,T,[],M)],'-',TOPT*[1;1],[-5;5],'x-k');
	set(gca,'ylim',[-20,20]);

	%semilogx(exp(X),[Y,polyval(P,X,[],M),polyval(P1,X,[],M)],'-',exp(T),[polyval(P,T,[],M),polyval(P1,T,[],M)],'-',TOPT*[1;1],[-5;5],'x-k');
end


return

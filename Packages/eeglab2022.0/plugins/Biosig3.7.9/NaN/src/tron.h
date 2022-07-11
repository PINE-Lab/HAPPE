/* 
  This code was extracted from liblinear-2.2.1 in Feb 2019 and 
  modified for the use with Octave and Matlab


Copyright (c) 2007-2019 The LIBLINEAR Project.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither name of copyright holders nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/ 

#ifndef _TRON_H
#define _TRON_H

class function
{
public:
	virtual double fun(double *w) = 0 ;
	virtual void grad(double *w, double *g) = 0 ;
	virtual void Hv(double *s, double *Hs) = 0 ;

	virtual int get_nr_variable(void) = 0 ;
	virtual void get_diag_preconditioner(double *M) = 0 ;
	virtual ~function(void){}
};

class TRON
{
public:
	TRON(const function *fun_obj, double eps = 0.1, double eps_cg = 0.1, int max_iter = 1000);
	~TRON();

	void tron(double *w);
	void set_print_string(void (*i_print) (const char *buf));

private:
	int trpcg(double delta, double *g, double *M, double *s, double *r, bool *reach_boundary);
	double norm_inf(int n, double *x);

	double eps;
	double eps_cg;
	int max_iter;
	function *fun_obj;
	void info(const char *fmt,...);
	void (*tron_print_string)(const char *buf);
};
#endif

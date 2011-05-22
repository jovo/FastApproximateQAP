function [f0new, salpha] = lines(type, x,d,g,A,B)
%function [f0new, salpha] = lines(type, x,d,g,A,B)
% line search: 0=> salpha=1, 1=> search
% 
% [f0new, salpha] = lines(type, T,O, x,d,n,m)
%--------------------
%n=length(dpi); m=length(dB(1,:)); 
%en=ones(n,1); em=ones(m,1);
% Louis J. Podrazik circa 1996
%
% IDA Center for Computing Sciences
%  (c) 1996, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%

nT=max(size(g));
Debug=0;

dxg  = d'*g; 
if(dxg > 0) 
  fprintf(1,'warning: Nonimproving Direction, <d,g> = %g\n', dxg);
  %d=-d;
end

if(type == 0)     %fprintf(1, 'line 0\n');
  salpha=1;
elseif(type == 1) %fprintf(1, 'line 1\n'); 
  %salpha=1.0;
%	searchtype = 1: bisection search
%-----------------------------------------------------------------------

	LARGE = 7.7E+77;
	l_TOL1  = 2.0E-12;
	l_iter  = 0;
	alpha       = 0.0;
	alpha_u     = 0.0;
	alpha_c     = 0.0;

	Done        = 0;
	nu          = 1.25;
	leps        = 0.2;

  	% alpha_c (CONSTRAINED SEARCH: Find Boundary)
  	% ==============================================================
     	alpha_c = LARGE ;
  	tdelta  = alpha_c;

	%---> Simplex constraint search --------------------------------
	  for j = 1:nT
	      if (d(j,1) < -l_TOL1)
	        tdelta = - x(j,1) / d(j,1);
	        alpha_c = min (tdelta, alpha_c);	% << COMPUTE ALPHA_C >>
	      end
	  end
	%---------------------------------------------------------------

	if(alpha_c == 0.0)
           fprintf(1,'alpha_c Error: alpha_c = %f\n',alpha_c);
	elseif ( alpha_c == LARGE ) 
           fprintf(1,'Constrained Line Search Error \n');
	end

	[F0_0]=fun(x,A,B);
	xt = x + alpha_c * d;
	[F0_c]=fun(xt,A,B);

	% ----> bisection search -----------------------------

	Done = 0;
%        write (6,*) 'bisection last_score = ', F0_0
  	alpha_u  = alpha_c;
	F0_R     = F0_c;
	alpha_n  = 1.0;

	l1 = alpha_c / 2.;
  	xt = x + l1 * d;
	[F0_L]=fun(xt,A,B);
	half = 2.00;

	while(~Done)
%	F0_L

	 l_iter = l_iter + 1;
	 if(l_iter > 25)
	   fprintf(1,'Error: Too many line searches = %d, alpha=%f\n',l_iter,l1 );
	   keyboard;
	 end
%         write (6,*) 'bisection score = ', F0_L, ' alpha_u=', l1
%         write (6,*) 'last score = ', F0_0,
%     >  	    ' current score =', F0_L,' alpha_u=', l1
%         write (6,*) 'F0_u = ', F0_u, ' l_iter=', l_iter
%        write(6,*)'L=',F0_L, ' R=',F0_R,' alpha_u=', l1,' last=',F0_0
	 if (F0_L < F0_R) 
	   F0_R     = F0_L;
	   alpha_u  = l1;
	   l1       = l1 / half;
  	   xt       = x + l1 * d;
	   [F0_L]=fun(xt,A,B);
	 else
	   if (F0_R < F0_0) 
	     F0_u = F0_R;
	     Done = 1;
	   else
	       l1 = l1 / half;
  	       xt = x + l1 * d;
	       [F0_L]=fun(xt,A,B);
	   end
	 end
	 % keyboard
	end

%======= cleanup ================================================

  	% compute alpha 
  	% ------------- 
    	if(F0_c <= F0_u)
	  F0        = F0_c;
    	  salpha     = alpha_c;
	  alpha_n   = 1.0;
    	  exact     = 0;
  	else
	  F0        = F0_u;
    	  salpha     = alpha_u;
	  alpha_n   = alpha_u / alpha_c;
    	  exact     = 1;
	end
    	mu = alpha_c;

%	GOTO 778 ! -->> Don't take boundary if improving point
%	!Take boundary if improving point --> MOD
%    	if(F0_c < F0_0) 
%	  F0        = F0_c;
%    	  salpha     = alpha_c;
%	  alpha_n   = 1.0;
%    	  exact     = 0;
%	end

        if (alpha < 0.0) 
	  fprintf(1,'+ error: alpha < 0 :alpha= %f\n)', alpha);
	  keyboard;
	end

        %if ( (iPRINT) | (Debug) ) 
          if (Debug) 
            fprintf(1,'linesearch iterations = %d\n', l_iter);
	  end
	  %end

    	tmp = abs( (F0-F0_0)/F0);
    	if( (tmp > 1.0e-12) && (F0 > F0_0) ) 
          fprintf(1,'Nonmonotone Line Search Error: cost > last cost, cost = %f, last cost = ',F0,F0_0);
	  keyboard;
	end


%-------------------------------
elseif (type==2)
% assume the function is quadratic
% derivative at alpha=0:
b = g' * d;
% constant term at alpha=0
c = fun(x,A,B);
% get second order coeff
fun_vertex = fun(x+d,A,B);
a = fun_vertex - b - c;
if( abs(a)<eps) 
   %disp('function is linear'); 
   salpha=1;
else
  salpha =min(1,max(-b/(2*a),0));
end;
fun_alpha = fun(x+salpha*d,A,B);
% check quadratic function
qfun_alpha = a*salpha*salpha + b*salpha + c;
if ((abs(a)>=eps)&&(abs(fun_alpha-qfun_alpha)>1000*abs(fun_alpha)*eps)) 
	fprintf(1,'quadratic search error: %g, %g\n',fun_alpha,qfun_alpha);
	keyboard;
end

if (fun_alpha>c) salpha=0; fun_alpha=c; end

if (fun_alpha>fun_vertex) salpha=1; fun_alpha = fun_vertex; end
f0new = fun_alpha;	

%-------------------------------

else
  fprintf(1, 'unsupported line\n'); keyboard;
end

xt=x + salpha*d;
f0new = fun(xt,A,B);

% salpha, f0new
% lmin= fmin('fun1dim',0,alpha_c,[0,1.e-6],x,d,T,O,n,m,scale)
% yy=fun1dim(lmin,x,d,T,O,n,m,scale)
% salpha=lmin; f0new=yy;

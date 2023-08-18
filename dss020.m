function ux=dss020(xl,xu,n,u,v)

%Routine dss020 implements fourth order directional differencing
%
%Intended for first order hyperbolic pde systems, e.g.
%
%       U +v*U  = 0
%        t    x
%
%Positive v for medium flowing in direction of higher x.
%
%Based on a routine with the same name in a book by W.E.Schiesser,
%The numerical method of lines, page 135.

dx    = (xu-xl)/(n-1);
rdx = 1/(12*dx);

if (v > 0)
        ux(n) = rdx*(  3*u(n-4)-16*u(n-3)+36*u(n-2)-48*u(n-1)+25*u(n));
        ux(1) = rdx*(-25*u(1)+48*u(2)-36*u(3)+16*u(4)- 3*u(5));
        ux(2) = rdx*(- 3*u(1)-10*u(2)+18*u(3)- 6*u(4)+ 1*u(5));
        ux(3) = rdx*(  1*u(1)- 8*u(2)+ 0*u(3)+ 8*u(4)- 1*u(5));
        for j=4:(n-1),
                ux(j) = rdx*(- 1*u(j-3)+ 6*u(j-2)-18*u(j-1)+10*u(j)+3*u(j+1));
        end
else

        ux(n)   = rdx*(  3*u(n-4)-16*u(n-3)+36*u(n-2)-48*u(n-1)+25*u(n));
        ux(n-1) = rdx*( -1*u(n-4)+ 6*u(n-3)-18*u(n-2)+10*u(n-1)+ 3*u(n));
        ux(n-2) = rdx*(  1*u(n-4)- 8*u(n-3)+ 0*u(n-2)+ 8*u(n-1)- 1*u(n));

        ux(1) = rdx*(-25*u(1)+48*u(2)-36*u(3)+16*u(4)- 3*u(5));
        for j=2:(n-3),
                %ux(j) = rdx*(- 1*u(j-1)+ 6*u(j)-18*u(j+1)+10*u(j+2)+3*u(j+3));
                ux(j) = rdx * (-3*u(j-1) -10*u(j) + 18*u(j+1) -6*u(j+2) + u(j+3) );
        end
end


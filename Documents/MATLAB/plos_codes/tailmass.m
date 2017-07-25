function check2 = tailmass( m, s, eps, T2, sd )
%TAILMASS Summary of this function goes here
%   Detailed explanation goes here



        % BEGIN FUNCTION TailMass
        % Input parameters: m, s, eps, T2, nu, sd,
        % Outputs: check2
            % to estimate the error (in calculating probabilities the of the
            % data), that results from approximating the first pdf as a
            % point-mass distribtion, find the maximum of the absolute value
            % of the derivative of the second pdf.
            gp=gp_max(m(2),s(2));

            % determine the radius, r, of a small interval over which the
            % 1. second pdf, g, is approximately constant, i.e. changes by less than eps/3
            % over any interval with that radius
            % and 2. g(t) is small for t<r.  ?????
            r=min(eps/(3*gp),T2);
            checkval=onestagepdf2(r,m(2),s(2));
            while checkval>=eps/2
                r=r/2;
                checkval=onestagepdf2(r,m(2),s(2));
            end

            % get the average value of the first pdf.  This is the point at
            % which its mass is concentrated.
            nu=1/m(1);

            % get the maximum value of the second pdf,f.
            gm=onestagepdf2(T2,m(2),s(2));

            % ????
            % Tu is the upper limit of integral for approximating
            % int_{r+nu}^{infty}g(s)ds. This is in latex.
            Tu=max(100,nu+r+1000*sd(1));

            % ????
            checkerror=100;
            hh=.001;
            numinusr=nu-r;
            nuplusr=nu+r;
            teeyou=Tu;
            %we are integrating the first pdf from zero to nu-r, and from nu+r to infinity to see if the tails
            %are small in probability.
            LeftTail=.001*sum(onestagepdf2((0:.001:numinusr),m(1),s(1)));
            RightTail=.001*sum(onestagepdf2((nuplusr:.001:teeyou),m(1),s(1)));
            check1=LeftTail+RightTail;
            %Reduce step size in above Riemann sum until the error is small,
            %meaning that the Riemann sum is converging.
            while checkerror>10^(-4)
                hh=.5*hh;
                ck1=hh*sum(onestagepdf2((0:hh:nu-r),m(1),s(1)))+hh*sum(onestagepdf2((nu+r:hh:Tu),m(1),s(1)));
                checkerror=abs(check1-ck1);
                check1=ck1;
            end
            check2=gm*check1;
            %     end
        % END FUNCTION TailMass



end


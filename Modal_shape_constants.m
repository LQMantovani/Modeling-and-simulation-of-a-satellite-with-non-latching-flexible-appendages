%Function to find lambdas e sigmas for the beam model used. This function considers the cantilevered-free beam with tip mass

global matrix Bodies nflex

for j=2:length(nflex) %To evaluate all the bodies
    if nflex(j)~=0
        L=Bodies.B(j).Length;
        r_e=Bodies.B(j).External_Radius;
        r_i=Bodies.B(j).Internal_Radius;
        rho=Bodies.B(j).Density;
        m=rho*L*pi*(r_e^2-r_i^2);
        Mt=Bodies.B(j).Tip_mass;
        Bodies.B(j).mass=m+Mt;
        It=Bodies.B(j).Tip_Inertia;
        
        %Lambda is determined from the detm=0. Since mathematica found no
        %analytical solution and it depends on the two bodies proprieties, it
        %is solved in the beggining of each simulation and stored in Bodies
        %structure
        detm=@(lm) (-1).*cos(lm).^2+(-1).*It.*L.^(-4).*lm.^4.*m.^(-2).*Mt.*cos( ...
            lm).^2+(-2).*cos(lm).*cosh(lm)+2.*It.*L.^(-4).*lm.^4.*m.^( ...
            -2).*Mt.*cos(lm).*cosh(lm)+(-1).*cosh(lm).^2+(-1).*It.*L.^( ...
            -4).*lm.^4.*m.^(-2).*Mt.*cosh(lm).^2+2.*It.*L.^(-3).*lm.^3.* ...
            m.^(-1).*cosh(lm).*sin(lm)+2.*L.^(-1).*lm.*m.^(-1).*Mt.* ...
            cosh(lm).*sin(lm)+(-1).*sin(lm).^2+(-1).*It.*L.^(-4).* ...
            lm.^4.*m.^(-2).*Mt.*sin(lm).^2+2.*It.*L.^(-3).*lm.^3.*m.^( ...
            -1).*cos(lm).*sinh(lm)+(-2).*L.^(-1).*lm.*m.^(-1).*Mt.*cos( ...
            lm).*sinh(lm)+sinh(lm).^2+It.*L.^(-4).*lm.^4.*m.^(-2).*Mt.* ...
            sinh(lm).^2;
        
        %Initializing the vectors
        lambda=zeros(1,Bodies.B(j).nflex); sigma=lambda; AC=zeros(size(lambda,2),2);
        
        i=1;
        c0=.1;
        while i<=size(lambda,2)
            lambda(i)=fsolve(detm,c0,optimset('Disp','off')); %Solving lambda
            if i==1
                if lambda(i)~=0
                    i=i+1; %If lambda is solved
                end
            else
                if lambda(i)>lambda(i-1)*1.1
                    i=i+1;
                end
            end
            c0=c0+.1; %Increase the intial guess to the next solution be found
        end
        for i=1:size(lambda,2)
            
            %Solving for sigma
            sigma(i)=(sin(lambda(i))-sinh(lambda(i))+lambda(i)*(Mt/(m*L))*(cos(lambda(i))-cosh(lambda(i))))...
                /(cos(lambda(i))+cosh(lambda(i))-lambda(i)*(Mt/(m*L))*(sin(lambda(i))-sinh(lambda(i))));
            
            %Solving for the AC constants
            lm=lambda(i);
            matrix=[cos(lm)+cosh(lm)+(-1).*It.*L.^(-3).*lm.^3.*m.^(-1).*(sin( ...
                lm)+sinh(lm)),It.*L.^(-3).*lm.^3.*m.^(-1).*(cos(lm)+(-1).* ...
                cosh(lm))+sin(lm)+sinh(lm);L.^(-1).*lm.*m.^(-1).*Mt.*(cos( ...
                lm)+(-1).*cosh(lm))+sin(lm)+(-1).*sinh(lm),(-1).*cos(lm)+( ...
                -1).*cosh(lm)+L.^(-1).*lm.*m.^(-1).*Mt.*(sin(lm)+(-1).*sinh( ...
                lm))];
            AC(i,:)=fsolve(@Modal_shapes_AC,[1;1],optimset('Disp','off'));
        end
        
        Bodies.B(j).lambda=lambda;
        Bodies.B(j).sigma=sigma;
        Bodies.B(j).AC=AC;
    end
end
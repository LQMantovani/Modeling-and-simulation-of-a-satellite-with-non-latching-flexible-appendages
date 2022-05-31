%FILE TO GENERATE SYSTEM GRAPHS

%Plot general configurations
size_l=25; %Font size
wd=2;      %Line width

nf=2; %Specify the figure number

%Main body position in the ECI
Xr = X(:,Bodies.B(1).pos_R);
raio0 = norm(Xr(1,1:3));
figure(nf)
plot3(Xr(:,1)/raio0, Xr(:,2)/raio0, Xr(:,3)/raio0); xlabel('X'); ylabel('Y'); zlabel('Z'); grid minor; %xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
nf = nf + 1;

%Bodies angular velocity
clear lgC;
figure(nf); vec_legend=['p^{%-.0f}' 'q^{%-.0f}' 'r^{%-.0f}'];
for j=1:3
    for i=1:length(nflex)
        subplot(3,1,j); plot(t,X(:,Bodies.B(i).pos_omega(j))*180/pi); xlabel(' t (s) '); ylabel(' degrees/s '); lgC{i,j}=num2str(i,vec_legend(1+9*(j-1):9*(j-1)+9));
        hold on; 
    end
    hold off;
    legend(lgC{:,j}); grid minor;
    sgtitle('Bodies angular velocity - BRF');
end
nf = nf + 1;

%Bodies angular position
if Bodies.System.quat == 0
    clear lgC vec_legend;
    figure(nf); vec_legend=['phi^{%-.0f}  ' 'theta^{%-.0f}' 'psi^{%-.0f}  '];
    for j=1:3
        parfor i=1:length(nflex)
            subplot(3,1,j); plot(t,X(:,Bodies.B(i).pos_Theta(j))*180/pi); xlabel(' t (s) '); ylabel(' degrees '); lgC{i,j}=num2str(i,vec_legend(1+13*(j-1):13*(j-1)+13));
            hold on;
        end
        hold off;
        legend(lgC{:,j}); grid minor;
        sgtitle('Bodies angular position - LVLH');
    end
    nf = nf + 1;
    
    clear lgC vec_legend
    figure(nf);
    subplot(311); plot(t,X(:,Bodies.B(1).pos_Theta(1))*180/pi); xlabel(' t (s) '); ylabel(' \phi (degrees)'); grid minor;
    subplot(312); plot(t,X(:,Bodies.B(1).pos_Theta(2))*180/pi); xlabel(' t (s) '); ylabel(' \theta (degrees)'); grid minor;
    subplot(313); plot(t,X(:,Bodies.B(1).pos_Theta(3))*180/pi); xlabel(' t (s) '); ylabel(' \psi (degrees)'); grid minor;
    sgtitle([Bodies.B(1).Name,' angular position (LVLH)']);
    nf = nf + 1;
else
    
    clear lgC vec_legend
    figure(nf);
    Xq = plot_degree(X(:,Bodies.B(1).pos_Theta));
    subplot(311); plot(t,Xq(:,1,1)*180/pi); xlabel(' t (s) '); ylabel('\phi'); grid minor;
    subplot(312); plot(t,Xq(:,2,1)*180/pi); xlabel(' t (s) '); ylabel('\theta'); grid minor;
    subplot(313); plot(t,Xq(:,3,1)*180/pi); xlabel(' t (s) '); ylabel('\psi'); grid minor;
    sgtitle([Bodies.B(1).Name,' angular position (LVLH)']);
    nf = nf + 1;
end

%Bodies elastic modes
clear lgC vec_legend;
vec_legend = ['qf^{%-.0f}_{%-.0f}']; vec_legend_vel=['dqf^{%-.0f}_{%-.0f}/dt'];
for j = 2:length(nflex)
    if nflex(j) ~= 0
        figure(nf);
        if(max(Bodies.B(j).Flexible_modes) < Bodies.B(j).nflex)
            n_flexibles = Bodies.B(j).Flexible_modes(Bodies.B(j).Flexible_modes > 0);
            Axis_labels = ['X^{%-.0f} axis'; 'Y^{%-.0f} axis'; 'Z^{%-.0f} axis'];
            Axis_labels = Axis_labels(Bodies.B(j).Flexible_modes > 0,:);
            
            for i = 1:n_flexibles(1)
                subplot(221); plot(t,X(:,Bodies.B(j).pos_qf(i))); xlabel(' t (s) '); ylabel(' qf '); lgC{i}=num2str([j,i],vec_legend); hold on;
            end
            title(num2str(j,Axis_labels(1,:)));
            legend(lgC); hold off; grid minor; clear lgC
            
            for i = 1:n_flexibles(2)
                subplot(222); plot(t,X(:,Bodies.B(j).pos_qf(i + n_flexibles(1)))); xlabel(' t (s) '); ylabel(' qf '); lgC{i}=num2str([j,i],vec_legend); hold on;
            end
            title(num2str(j,Axis_labels(2,:)));
            legend(lgC); hold off; grid minor; clear lgC
            
            for i = 1:n_flexibles(1)
                subplot(223); plot(t,X(:,Bodies.B(j).pos_qfp(i))); xlabel(' t (s) '); ylabel(' dqf/dt '); lgC_v{i}=num2str([j,i],vec_legend_vel); hold on;
            end
            legend(lgC_v); hold off; grid minor; clear lgC_v
            
            for i = 1:n_flexibles(2)
                subplot(224); plot(t,X(:,Bodies.B(j).pos_qfp(i + n_flexibles(1)))); xlabel(' t (s) '); ylabel(' dqf/dt '); lgC_v{i}=num2str([j,i],vec_legend_vel); hold on;
            end
            legend(lgC_v); hold off; grid minor; clear lgC_v
            
            sgtitle('Bodies flexible modes');
        else
            for i=1:nflex(j)
                subplot(211); plot(t,X(:,Bodies.B(j).pos_qf(i))); xlabel(' t (s) '); ylabel(' qf '); lgC{i}=num2str([j,i],vec_legend); hold on;
            end
            legend(lgC); hold off; grid minor;
            for i=1:nflex(j)
                subplot(212); plot(t,X(:,Bodies.B(j).pos_qfp(i))); xlabel(' t (s) '); ylabel(' dqf/dt '); lgC_v{i}=num2str([j,i],vec_legend_vel); hold on;
            end
            legend(lgC_v); hold off; grid minor;
            sgtitle([Bodies.B(j).Name, ' flexible modes']);
        end
        nf = nf + 1;
    end
end


plot_relative_data
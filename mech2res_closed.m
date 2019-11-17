function [serpentine_res,stretch_array] = mech2res_closed(num_beams,prestretch,stretch_int)
%inputs num_beams = half the number of beams (symm), num_segs, prestrain,disp_int,tar_disp

stretch_array = [1:stretch_int:1.02*prestretch];

%%%%% Beam geometry properties %%%%%
a = 10e-6; %Length of short part of meander
w = 200e-6-a; %beam width
b = 10e-3-w; %Length of long part of meander

prestretch_fcn = prestretch;

[EI, Rs, Rl, Rc, q_eff]= beam_properties(a,b,w);

serpentine_res = zeros(1,length(stretch_array));

[prestretch_pos,beam_ledge] = prestretch_calc(prestretch_fcn,num_beams,a,w);
    
    for stretch_count = 1:length(stretch_array)
        %figure;
        curr_stretch = stretch_array(stretch_count); %Get the stretch applied to the sensor from the stretch array
        
        eq_pos = prestretch_pos*curr_stretch; %Recalculate the equilibrium position (zero energy stored in elastic membrane) 
        
        [~,L_eff] = find_Umin(EI,q_eff,w,num_beams,b,a,eq_pos,beam_ledge); %Find the beam position that results in minimum stored energy
                
        serpentine_res(stretch_count) = area_to_resistance(num_beams,b,Rs,Rl,Rc,L_eff); %Calculate resistance of serpentine based on beam positions
    end

figure;
plot(100*(stretch_array-1),serpentine_res,'.') %Plot resistance vs. stretch
stretch_array = 100.*(stretch_array-1); %Change to match stretch array values with collected data
end

function [EI, Rs, Rl, Rc, q_eff]= beam_properties(a,b,w)
%Beam properties and dimensions
t = 145e-6; %beam thickness
E = 200e9; %Young's modulus
EI = E*(t*w^3)./12;

k_eff = 5e3; %Spring constant of the elastic membrane
q_eff = k_eff/b; %Distributed spring constant

%%%%% Electrical properties %%%%%

%Values for anisotopy vary by structure of carbon fiber layup. Uncomment
%desired values

%0-90-0-90-0 pitch
% rho_90 = 89.427e-6;
% rho_0 = rho_90/1.24;

% %0-90-0-90-0 width
%  rho_90 = 73.421e-6;
%  rho_0 = rho_90/1.24;

% Anisotropy low
%rho_90 = 88.050e-6;
%rho_0 = rho_90/3.8;

% Anisotropy high
rho_90 = 273.68e-6;
rho_0 = rho_90/3.8;

Rc = 1e-1; %Contact resistance 
Rs = rho_0*(a+w)/(w*t); %Calculation for resistance of short component of serpentine
Rl = rho_90*(b+w)/(w*t); %Calculation for resistance of long component of serpentine
end

function [prestretch_pos,beam_ledge_base] = prestretch_calc(prestretch_fcn,num_beams,a,w)
beam_ledge_base = a/2+(a+w)*[0:1:num_beams-1];
prestretch_val = polyval(prestretch_fcn,beam_ledge_base);
prestretch_pos = beam_ledge_base./prestretch_val;
end


function [beam_ledge_base,L_eff] = find_Umin(EI,q_eff,w,num_beams,b,a,eq_pos,beam_ledge_base)

L_eff= b*ones(1,num_beams); %Assume the effective length is the beam length (no contact)

beam_ledge_tip = beam_ledge_base; %Set the beam curvature to 0

%The code can be modified to plot the position of each beam as a fill. This code is normally commented out to speed up execution 
% truss_x_up = b*ones(1,5);
% truss_x_down = zeros(1,5);
% x_pos = [0:100e-6:b];
% x_plot_odd = [x_pos truss_x_up b-x_pos truss_x_down];
% x_plot_even = [b-x_pos truss_x_down x_pos truss_x_up];
% 
% D = 0; %D represents the global location of the base of the beam
% 
% Leff_prev = b*ones(1,num_beams);

for beam_num = 1:num_beams
    
    %Find displacement of current beam that minimizes substrate energy
    subs_disp = sum(eq_pos(beam_num:end)-beam_ledge_base(beam_num:end))/(num_beams+1-beam_num);
        
    %Find displacement to left (closer to symmetry) neighboring beam
    if beam_num > 1
        gapLR = beam_ledge_base(beam_num-1)+w-beam_ledge_tip(beam_num);
    else
        gapLR = -a/2;
    end
        
    if gapLR > subs_disp 
        %If |displacement| to previous beam is smaller than displacement to eq_pos AND disp to eq_pos is negative
    
        beam_disp = gapLR; %Beam will move until it contacts next beam (distance of gapLR)
        
        %dU_vec_point = [q_eff/2*(eq_pos(beam_num)^2-2*eq_pos(beam_num)*(D+beam_disp/2)+(beam_disp^2*13/35)+beam_disp*D+D^2),0,0,0,-18*EI*beam_disp^2];
        dU_vec_dist = [q_eff/2*(eq_pos(beam_num)^2-2*eq_pos(beam_num)*(D+4*beam_disp/15)+((beam_disp/3)^2*88/63)+8*D*beam_disp/15+D^2),0,0,0,-16*11*EI*beam_disp^2/5];

        L_eff_sol = roots(dU_vec_dist);

        for sol_num=1:length(L_eff_sol)
            if isreal(L_eff_sol(sol_num)) 
                if L_eff_sol(sol_num) > 0 
                    if L_eff_sol(sol_num) < b
                        L_eff(beam_num) = L_eff_sol(sol_num);
                    end
                end
            end
        end
        
        

        tip_disp = beam_disp;
         
    else
        %If the gap displacement is more negative OR gap displacement is
        %negative and beam displacement for equilibrium position is positive... 
        tip_disp = subs_disp;   
    end
        
    %Move tip end of beam by whatever the tip displacement is   
    beam_ledge_tip(beam_num) = beam_ledge_base(beam_num) + tip_disp;    
        
    %Update position of following beams
    if beam_num < num_beams
        beam_ledge_base(beam_num+1:num_beams) = beam_ledge_base(beam_num+1:num_beams) + tip_disp;
        beam_ledge_tip(beam_num+1:num_beams) = beam_ledge_tip(beam_num+1:num_beams) + tip_disp;
    end
    
   
    
       %%%%% Plotting fill here %%%%%
%     x_up = [0:100e-6:L_eff(beam_num)];
%     
%     %length(x_up)
%     y_up = ones(1,ceil(b./100e-6))*tip_disp;
%     y_up(1:length(x_up)) = tip_disp/(3*L_eff(beam_num)^4)*(4*L_eff(beam_num)*x_up.^3-x_up.^4);
%     
%     y_down = max(y_up) + w -flipud(y_up);
%     
%     truss_y_up = max(y_up) + [0:w/4:w];
%     truss_y_down = max(y_up) +w-[0:w/4:w];
%     y_plot = [y_up truss_y_up y_down truss_y_down] + beam_ledge_base(beam_num);
    
%     hold on;
%     if mod(beam_num,2) == 1 
%         fill(y_plot,x_plot_odd,'k')
%         plot(y_plot,x_plot_odd,'r')
%     else
%         fill(y_plot,x_plot_even,'k')
%         plot(y_plot,x_plot_even,'r')
%     end
    %%%%%

    
end
% L_eff
% hold on;
% plot(L_eff-Leff_prev);
% L_eff = Leff_prev;
end




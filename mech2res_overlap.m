function [serpentine_res,strain_array] = mech2res_overlap(num_beams,num_segs, prestrain,strain_int)
%inputs num_beams = half the number of beams (symm), num_segs, prestrain,disp_int,tar_disp

tic
strain_array = [0:strain_int:1.3*prestrain];

[EI, a, b, w, Rs, Rl, rho_c, q_eff, xi_int, xi_array]= beam_properties(num_segs);

[beam_ledge,norm_disp] = make_disp_arrays(a,w,EI,xi_array,num_beams);
prestrain_pos = beam_ledge;

serpentine_res = zeros(1,length(strain_array));

for strain_count = 1:length(strain_array)
    currstrain = strain_array(strain_count);
    
    %Find position of all beams once displacement equilibrium is reached
    eq_pos = prestrain_pos.*((100+currstrain)/(100+prestrain)); %Position of zero energy stored in substrate
    
    [beam_ledge,overlap] = find_Umin_wstrain(EI,q_eff,xi_int,w,num_beams,xi_array,eq_pos,prestrain_pos,norm_disp);
    %beam_plot(beam_ledge,strain_count,strain_int,b,xi_array,num_beams,w,a);
    contact_res = calc_overlap_res(rho_c, overlap,num_beams);
    serpentine_res(strain_count) = area_to_resistance_w_contact(num_beams,Rs,Rl,contact_res);
end
toc
plot(strain_array,serpentine_res,'.')
end


function [EI, a, b, w, Rs, Rl, rho_c, q_eff, xi_int, xi_array]= beam_properties(num_segs)
%Beam properties and dimensions
t = 140e-6; %beam thickness
w = 195e-6; %beam width
E = 200e9; %Young's modulus
EI = E*(t*w^3)./12;
a = 5e-6; %Length of short part of meander
b = 19.8e-3; %Length of long part of meander

k_eff = (2500); %Spring constant of one gap section of the substrate (need to measure)
q_eff = k_eff/b;

%%%%% Beam setup variables %%%%%
%num_segs = 9; %# of segments in each beam
xi_int = b/num_segs; %Length of ea  ch segment
xi_array = [0:xi_int:b]; %Location of each segment (base is 0, tip is b)

%%%%% Electrical properties %%%%%
rho_90 = 7.17e-5;
rho_0 = rho_90/4;
rho_c = 1e-2;
Rs = rho_0*a/(w*t);
Rl = rho_90*b/(w*t);

end

function [beam_ledge,norm_disp] = make_disp_arrays(a,w,EI,xi_array,num_beams)
%%%%%  Initialize beam positions %%%%%

beam_ledge = [a/2+[0:1:num_beams-1]*(a+w)]'*ones(1,length(xi_array)); %Initial, prestrained beam position

norm_disp_interval = 0.05;
norm_disp = ones(length(xi_array)+floor(1/norm_disp_interval)-1,length(xi_array));

%Add "saturated" beam displacements to the norm_disp array
for xi_ind = 1:length(xi_array)-1
    base_xi = xi_array(1:xi_ind+1);
    norm_disp(xi_ind,1:length(base_xi)) = -(base_xi).^3./(6*EI)+base_xi(xi_ind+1).*((base_xi).^2)./(4*EI);
    norm_disp(xi_ind,1:length(base_xi)) = norm_disp(xi_ind,1:length(base_xi))/max(norm_disp(xi_ind,1:length(base_xi)));
end

%Add "unsaturated" beam displacements to the norm_disp array
norm_disp(length(xi_array):length(xi_array)+floor(1/norm_disp_interval),:) = [1:-norm_disp_interval:0]'*norm_disp(length(xi_array)-1,:);
norm_disp(length(xi_array)+floor(1/norm_disp_interval)+1:length(xi_array)+2*floor(1/norm_disp_interval),:) = -1*[norm_disp_interval:norm_disp_interval:1]'*norm_disp(length(xi_array)-1,:);
end

function [beam_ledge,overlap] = find_Umin_wstrain(EI,q_eff,xi_int,w,num_beams,xi_array,eq_pos,beam_ledge,norm_disp)

overlap = zeros(num_beams,length(xi_array));

for beam_num = 1:num_beams
    
    %%%%%  Calculate the initial moment and displacement of each beam by finding minimum energy  %%%%%
    if beam_num == 1
        beam_redge = zeros(1,length(xi_array)); %Right edge is at symm for 1st beam
    else
        beam_redge = fliplr(beam_ledge(beam_num-1,:))+w; %Right edge is w+prev beam left edge. Flip to match base/tip
    end
    
    max_gap = -1*max(beam_ledge(beam_num,:)-beam_redge);
    
    guess_disp = norm_disp.*max_gap; %Potential displacement array 
    
    %Calculate initial position from energy minimization
    angle = gradient(bsxfun(@plus,guess_disp,beam_ledge(beam_num,:)),xi_int); %Find y' given displacement
    M_beam = gradient(angle,xi_int).*EI; %Find moment given y'
    
    U_strain = 1/(2*EI)*max(cumtrapz(xi_array,(M_beam.^2)'))'; %Find beam strain  
    
    U_subs = 0.5*q_eff*xi_int*sum((bsxfun(@minus,eq_pos(beam_num,:)-beam_ledge(beam_num,:),guess_disp)).^2,2); %Energy stored in substrate-beam as fcn of guess_disp
    
    for beams = beam_num+1:num_beams
        beam_move = guess_disp(:,end)*ones(1,length(xi_array));
        U_subs = U_subs + 0.5*q_eff*xi_int*sum((bsxfun(@minus,eq_pos(beams,:)-beam_ledge(beams,:),beam_move)).^2,2); %Energy stored in substrate-beam as fcn of guess_disp
    end
    
    [~,U_min_index]=min(U_strain+U_subs); %Find guess_disp that yields min energy
    
    %%%Update beam position based on which beam position gives min energy
    new_beam_ledge = max(beam_ledge(beam_num,:)+guess_disp(U_min_index,:),beam_redge);
    overlap(beam_num,:) = -1*min(beam_ledge(beam_num,:)+guess_disp(U_min_index,:)-beam_redge,0);

    delta_disp = new_beam_ledge(end) - beam_ledge(beam_num,end); %Actual change in displacement with beam contact condition applied
    
    beam_ledge(beam_num,:) = new_beam_ledge;
    
    %Update subsequent beam positions
    if beam_num < num_beams
        beam_ledge(beam_num+1:num_beams,:) = beam_ledge(beam_num+1:num_beams,:)+delta_disp;
    end    
end
end

function beam_plot(beam_ledge,disp_count,disp_int,b, xi_array,num_beams,w,a)
figure;
hold on;
title(strcat('Displacement = ',num2str((disp_count)*disp_int)));
xlabel('X position (m)')
ylabel('Y position (m)')
for beam_num = 1:num_beams
    if  mod(beam_num,2)==1
        y_data = xi_array;
        truss_y_loc = [b b];
    else
        y_data = b-xi_array;
        truss_y_loc = [0 0];
    end
    plot(beam_ledge(beam_num,:),y_data,'k')
    plot(beam_ledge(beam_num,:)+w,y_data,'g')
    plot([beam_ledge(beam_num,end)+w beam_ledge(beam_num,end)+a+w],truss_y_loc,'k')
end
end

function contact_res = calc_overlap_res(rho_c, overlap,num_beams)
contact_res = zeros(1,num_beams);
for beam_num = 1:num_beams
    cond_c = sum(1./(rho_c./overlap(beam_num,:)+100));
    contact_res(beam_num) = 1/cond_c;
    if contact_res(beam_num) == -Inf
        contact_res(beam_num) = 1e9;
    end
end
end
    
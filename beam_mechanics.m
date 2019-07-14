function beam_mechanics(num_beams,num_segs, prestrain,disp_int,tar_disp)
%inputs num_beams = half the number of beams (symm), num_segs

%num_beams = 5;
%num_segs = 10;
%prestrain = 10; %Prestrain in %
%disp_int = 2.5e-6;
num_disp_loops = ceil(tar_disp/disp_int); %Target final displacement

[EI, a, b, w, q_eff, xi_int, xi_array]= beam_properties(num_segs);

[beam_ledge,eq_pos,norm_disp] = make_disp_arrays(a,w,prestrain,EI,xi_array,num_beams);

for disp_count = 0:num_disp_loops
    
    %Find initial position of all beams once prestrain is released
    if disp_count == 0
        [beam_ledge,~] = find_Umin_pos_init(EI,q_eff,xi_int,w,num_beams,xi_array,eq_pos,beam_ledge,norm_disp);
        beam_plot(beam_ledge,disp_count,disp_int,b, xi_array,num_beams,w,a);
    
    %displace each beam by disp_int, calculate Umin and permanently displace lowest energy condition    
    else
        U_tot = zeros(1,num_beams);
        temp_beam_ledge_array = zeros(num_beams^2,length(xi_array));
        temp_eq_pos_array = zeros(num_beams^2,length(xi_array));
        
        for disp_beam_num = 1:num_beams
            %initialize temporary displacement arrays
            temp_beam_ledge = beam_ledge;
            temp_eq_pos = eq_pos;
            
            %displace beams from disp_beam_num onward by disp_int
            forced_beam_disp = norm_disp(length(xi_array)-1,:).*disp_int;
            
            temp_beam_ledge(disp_beam_num,:) = temp_beam_ledge(disp_beam_num,:) + forced_beam_disp;

            for beams = disp_beam_num+1:num_beams
                temp_beam_ledge(beams,:) = temp_beam_ledge(beams,:) + forced_beam_disp*mod(beams+1-disp_beam_num,2) + fliplr(forced_beam_disp)*mod(beams-disp_beam_num,2);
                temp_eq_pos(beams,:) = temp_eq_pos(beams,:)         + forced_beam_disp*mod(beams+1-disp_beam_num,2) + fliplr(forced_beam_disp)*mod(beams-disp_beam_num,2);
            end
            
            %find displacement of beams following beam that is moved
            [temp_beam_ledge,U_out] = Umin_w_disp(EI,q_eff,xi_int,w,num_beams,xi_array,temp_eq_pos,temp_beam_ledge,norm_disp,disp_beam_num);
            
            %add displacement/eq_pos to array
            temp_beam_ledge_array((disp_beam_num-1)*num_beams+1:(disp_beam_num)*num_beams,:) = temp_beam_ledge;
            temp_eq_pos_array((disp_beam_num-1)*num_beams+1:(disp_beam_num)*num_beams,:) = temp_eq_pos;
            
            %add minimum energy to array
            U_tot(disp_beam_num) = U_out;
        end
        
        %find the displacement of the beam that causes minimum energy
        [~,pos_min_index] = min(U_tot);
        
        %update beam positions
        beam_ledge = temp_beam_ledge_array((pos_min_index-1)*num_beams+1:(pos_min_index)*num_beams,:);
        
        %update equilibrium positions
        eq_pos = temp_eq_pos_array((pos_min_index-1)*num_beams+1:(pos_min_index)*num_beams,:);
        
        %plot updated beam positions
        if mod(disp_count,10) == 1
        beam_plot(beam_ledge,disp_count,disp_int,b, xi_array,num_beams,w,a);
        end
    end
end
end


function [EI, a, b, w, q_eff, xi_int, xi_array]= beam_properties(num_segs)
%Beam properties and dimensions
t = 90e-6; %beam thickness
w = 100e-6; %beam width
E = 100e9; %Young's modulus
EI = E*(t*w^3)./12;
a = 5e-6; %Length of short part of meander
b = 20e-3; %Length of long part of meander
%%%%%delam = 0;

%num_beams = 7; %Half number of beams.
k_eff = (2500); %Spring constant of one gap section of the substrate (need to measure)
q_eff = k_eff/b;

%%%%% Beam setup variables %%%%%
%num_segs = 9; %# of segments in each beam
xi_int = b/num_segs; %Length of ea  ch segment
xi_array = [0:xi_int:b]; %Location of each segment (base is 0, tip is b)

end

function [beam_ledge,eq_pos,norm_disp] = make_disp_arrays(a,w,prestrain,EI,xi_array,num_beams)
%%%%%  Initialize beam positions %%%%%

beam_ledge = [a/2+[0:1:num_beams-1]*(a+w)]'*ones(1,length(xi_array)); %Initial, prestrained beam position
eq_pos = (beam_ledge./(1+prestrain/100)); %Position of zero energy stored in substrate

norm_disp_interval = 0.1;
norm_disp = ones(length(xi_array)+floor(1/norm_disp_interval)-1,length(xi_array));

%Add "saturated" beam displacements to the norm_disp array
for xi_ind = 1:length(xi_array)-1
    base_xi = xi_array(1:xi_ind+1);
    norm_disp(xi_ind,1:length(base_xi)) = -(base_xi).^3./(6*EI)+base_xi(xi_ind+1).*((base_xi).^2)./(4*EI);
    norm_disp(xi_ind,1:length(base_xi)) = norm_disp(xi_ind,1:length(base_xi))/max(norm_disp(xi_ind,1:length(base_xi)));
end

%Add "unsaturated" beam displacements to the norm_disp array
norm_disp(length(xi_array):length(xi_array)+floor(1/norm_disp_interval)-1,:) = [1-norm_disp_interval:-1*norm_disp_interval:0]'*norm_disp(length(xi_array)-1,:);
end

function [beam_ledge,U_min] = find_Umin_pos_init(EI,q_eff,xi_int,w,num_beams,xi_array,eq_pos,beam_ledge,norm_disp)

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
    
    [U_min,U_min_index]=min(U_strain+U_subs); %Find guess_disp that yields min energy
    
    %%%Update beam position based on which beam position gives min energy
    new_beam_ledge = max(beam_ledge(beam_num,:)+guess_disp(U_min_index,:),beam_redge);
    
    delta_disp = new_beam_ledge(end) - beam_ledge(beam_num,end); %Actual change in displacement with beam contact condition applied
    
    beam_ledge(beam_num,:) = new_beam_ledge;
    
    %Update subsequent beam positions
    if beam_num < num_beams
        beam_ledge(beam_num+1:num_beams,:) = beam_ledge(beam_num+1:num_beams,:)+delta_disp;
    end    
end
end

function [beam_ledge,U_tot] = Umin_w_disp(EI,q_eff,xi_int,w,num_beams,xi_array,eq_pos,beam_ledge,norm_disp,start_beam)

U_tot = 0;

for beam_num=1:start_beam
    %Find the energy of each "fixed" beam (beams before displaced beam)
    U_tot = U_tot + find_Umin_prev_beam(EI,q_eff,xi_int,xi_array,eq_pos(beam_num,:),beam_ledge(beam_num,:));
end

for beam_num=start_beam+1:num_beams
    %Find the position and energy of each beam that can move to minimize energy
    [beam_ledge,U_min] = find_Umin_pos_succ(EI,q_eff,xi_int,w,num_beams,xi_array,eq_pos,beam_ledge,norm_disp,beam_num);
    U_tot = U_tot+U_min;
end
end

function U_tot = find_Umin_prev_beam(EI,q_eff,xi_int,xi_array,sing_eq_pos,sing_beam_ledge)

angle = gradient(sing_beam_ledge,xi_int); %Find angle given displacement
M_beam = gradient(angle,xi_int).*EI; %Find moment given angle
U_strain = 1/(2*EI)*max(cumtrapz(xi_array,(M_beam.^2)'))'; %Find beam strain

U_subs = sum(0.5*q_eff*xi_int*sum((sing_eq_pos-sing_beam_ledge).^2));
U_tot = U_strain+U_subs;
end

function [beam_ledge,U_min] = find_Umin_pos_succ(EI,q_eff,xi_int,w,num_beams,xi_array,eq_pos,beam_ledge,norm_disp,curr_beam)

%%%%%  Calculate the initial moment and displacement of each beam by finding minimum energy  %%%%%

beam_redge = fliplr(beam_ledge(curr_beam-1,:))+w; %Right edge is w+prev beam left edge. Flip to match base/tip

max_gap = -1*max(beam_ledge(curr_beam,:)-beam_redge);

guess_disp = norm_disp.*max_gap; %Potential displacement array

%Calculate initial position from energy minimization
angle = gradient(bsxfun(@plus,guess_disp,beam_ledge(curr_beam,:)),xi_int); %Find y' given displacement
M_beam = gradient(angle,xi_int).*EI; %Find moments given y'

U_strain = 1/(2*EI)*max(cumtrapz(xi_array,(M_beam.^2)'))'; %Find beam strain by integrating moment from base to tip

U_subs = 0.5*q_eff*xi_int*sum((bsxfun(@minus,eq_pos(curr_beam,:)-beam_ledge(curr_beam,:),guess_disp)).^2,2); %Energy stored in substrate-beam as fcn of guess_disp

U_subs_succ = 0;
for beams = curr_beam+1:num_beams
    beam_move = guess_disp*mod(beams+1-curr_beam,2)+fliplr(guess_disp)*mod(beams+1-curr_beam,2); %Need to fix beam move? HELP
    U_subs_succ = U_subs_succ + 0.5*q_eff*xi_int*sum((bsxfun(@minus,eq_pos(beams,:)-beam_ledge(beams,:),beam_move)).^2,2); %Energy stored in substrate-beam as fcn of guess_disp
end

[~,U_min_index]=min(U_strain+U_subs+U_subs_succ); %Find guess_disp that yields min energy for all subsequent beams

U_min = U_strain(U_min_index)+U_subs(U_min_index); %Return energy of just current beam/substrate

%%%Update beam position based on which beam position gives min energy
new_beam_ledge = max(beam_ledge(curr_beam,:) + guess_disp(U_min_index,:),beam_redge); %Check for overlap

delta_disp = new_beam_ledge-beam_ledge(curr_beam,:);

beam_ledge(curr_beam,:) = new_beam_ledge;

%Update subsequent beam positions
for beam_shift=curr_beam+1:num_beams
    %beam_ledge(beam_shift,:) = beam_ledge(beam_shift,:) +guess_disp(U_min_index,:)*mod(beam_shift+1-curr_beam,2) + fliplr(guess_disp(U_min_index,:))*mod(beam_shift-curr_beam,2); %Is correct with displacement?
    beam_ledge(beam_shift,:) = beam_ledge(beam_shift,:) + delta_disp*mod(beam_shift+1-curr_beam,2)+fliplr(delta_disp)*mod(beam_shift-curr_beam,2);
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

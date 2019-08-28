function [serpentine_res,stretch_array] = mech2res_closed(num_beams,prestretch,stretch_int)
%inputs num_beams = half the number of beams (symm), num_segs, prestrain,disp_int,tar_disp

stretch_array = [1:stretch_int:1.2*prestretch];

%%%%% Beam geometry properties %%%%%
a = 7.5e-6; %Length of short part of meander
w = 200e-6-a; %beam width
b = 10e-3-w; %Length of long part of meander

prestretch_fcn = prestretch;

[EI, Rs, Rl, Rc, q_eff]= beam_properties(a,b,w);

serpentine_res = zeros(1,length(stretch_array));

[prestretch_pos,beam_ledge] = prestretch_calc(prestretch_fcn,num_beams,a,w);

    for stretch_count = 1:length(stretch_array)
        curr_stretch = stretch_array(stretch_count);
        
        eq_pos = prestretch_pos*curr_stretch;
        
        [~,L_eff] = find_Umin(EI,q_eff,w,num_beams,b,a,eq_pos,beam_ledge);
        
        serpentine_res(stretch_count) = area_to_resistance(num_beams,b,Rs,Rl,Rc,L_eff);
    end

plot(100*(stretch_array-1),serpentine_res,'o')
end

function [EI, Rs, Rl, Rc, q_eff]= beam_properties(a,b,w)
%Beam properties and dimensions
t = 140e-6; %beam thickness
E = 200e9; %Young's modulus
EI = E*(t*w^3)./12;

E_TPU = 130e6;
t_TPU = 50e-6;

k_eff = 5000;%E_TPU*t_TPU; %Check these values with Seun
q_eff = k_eff/b;

%%%%% Electrical properties %%%%%
rho_90 = 75e-6;
rho_0 = rho_90/3.75;
Rc = 5*10e-3; %fitted parameter
Rs = rho_0*(a+w)/(w*t);
Rl = rho_90*(b+w)/(w*t);
end

function [prestretch_pos,beam_ledge] = prestretch_calc(prestretch_fcn,num_beams,a,w)
beam_ledge = a/2+(a+w)*[0:1:num_beams-1];
prestretch_val = polyval(prestretch_fcn,beam_ledge);
prestretch_pos = beam_ledge./prestretch_val;
end


function [beam_ledge_base,L_eff] = find_Umin(EI,q_eff,w,num_beams,b,a,eq_pos,beam_ledge)

L_eff= (b-1e-12)*ones(1,num_beams);

beam_ledge_tip = beam_ledge;
beam_ledge_base = beam_ledge;

for beam_num = 1:num_beams
    
    %Find displacement of current beam that minimizes substrate energy
    subs_disp = sum(eq_pos(beam_num:end)-beam_ledge_tip(beam_num:end))/(num_beams+1-beam_num);
    
    %Find displacement to left (closer to symmetry) neighboring beam
    
        gapLR = -a/2; 
    
    
    if gapLR > subs_disp 
        %If displacement to previous beam is |smaller| than displacement to eq_pos AND disp to eq_pos is negative
    
        beam_disp = gapLR;

        D = beam_ledge(beam_num);

        dU_vec_point = [q_eff/2*(eq_pos(beam_num)^2-2*eq_pos(beam_num)*(D+beam_disp/2)+(beam_disp^2*13/35)+beam_disp*D+D^2),0,0,0,-18*EI*beam_disp^2];
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

        %Update tip position to the beam displacement OR the right edge of
        %the previous beam, whichever is further away from symmetry
 
        if beam_num > 1
            tip_disp = max(beam_disp,beam_ledge_base(beam_num-1)+w-beam_ledge_tip(beam_num));
        else
            tip_disp = beam_disp;
        end
         
    else
        %If the gap displacement is more negative OR gap displacement is
        %negative and beam displacement for equilibrium position is positive... 
        tip_disp = subs_disp;   
    end
    
    beam_ledge_tip(beam_num) = beam_ledge_tip(beam_num) + tip_disp;
    
    %Update position of following beams
    if beam_num < num_beams
        beam_ledge_base(beam_num+1:num_beams) = beam_ledge_base(beam_num+1:num_beams) + tip_disp;
        beam_ledge_tip(beam_num+1:num_beams) = beam_ledge_base(beam_num+1:num_beams);
    end

end
end




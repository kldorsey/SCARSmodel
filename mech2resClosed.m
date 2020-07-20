function [serpentineRes,strainArray] = mech2resClosedVar(numBeams,prestretch,stretch_int,pitch,width,rho_0,rho_90,plotOn)
    %[serpentine_res,strain_array] = mech2res_closed_var(numBeams,prestretch,stretch_int,pitch,width,rho_0,rho_90,plotOn)

    %%%%% Effective beam geometry properties %%%%%
    gap = 10e-6; 
    pitch = pitch-gap; %change pitch to effective beam width
    width = width-pitch; %change width to effective beam length
    
    [EI, Rs, Rl, Rc, qEff]= beam_properties(gap,width,pitch,rho_0,rho_90);
    
    %%%%% Create array of stretches and find initial position
    stretchArray = [1:stretch_int:1.02*prestretch];
    serpentineRes = zeros(1,length(stretchArray));

    
    %Calculate beam pre-stretch position based on the prestretch function
    [prestretchPos,beamLedge] = prestretchCalc(prestretch,numBeams,gap,pitch);

    %For every stretch in the array, find the sensor resistance
    for stretchCount = 1:length(stretchArray) 

        %Get the stretch applied to the sensor from the stretch array
        currStretch = stretchArray(stretchCount); 
 
        %Recalculate the equilibrium position (zero energy stored in elastic membrane) 
        eqPos = prestretchPos*currStretch; 
        
        %Find the beam position that results in minimum stored energy
        [~,L_eff] = find_Umin(EI,qEff,pitch,numBeams,width,gap,eqPos,beamLedge); 
        
        %Calculate resistance of serpentine based on beam positions
        serpentineRes(stretchCount) = area_to_resistance(numBeams,width,Rs,Rl,Rc,L_eff); 
    end
    
    strainArray = 100.*(stretchArray-1);     %Strain array in percentage
    
    if plotOn == 1
        figure;
        plot(strainArray,serpentineRes,'.') %Plot resistance vs. stretch
    end
end

function [EI, Rs, Rl, Rc, q_eff]= beam_properties(a,width,pitch,rho_0,rho_90)
	%Beam properties and dimensions
	t = 145e-6; %beam thickness
	E = 200e9; %Young's modulus
	EI = E*(t*pitch^3)./12;

	k_eff = 5e3; %Spring constant of the elastic membrane
	q_eff = k_eff/width; %Distributed spring constant

	Rc = 1e-1; %Contact resistance 
	Rs = rho_0*(a+pitch)/(pitch*t); %Calculation for resistance of short component of serpentine
	Rl = rho_90*(width+pitch)/(pitch*t); %Calculation for resistance of long component of serpentine
end

function [prestretchPos,beamLedgeBase] = prestretchCalc(prestretch,numBeams,gap,pitch)
	beamLedgeBase = gap/2+(gap+pitch)*[0:1:numBeams-1];
	prestretchPos = beamLedgeBase./prestretch;
end


function [beam_ledge_base,L_eff] = find_Umin(EI,q_eff,pitch,num_beams,width,gap,eq_pos,beam_ledge_base)

	L_eff= width*ones(1,num_beams); %Assume the effective beam length is the sensor width (i.e., no contact)
	beam_ledge_tip = beam_ledge_base; %Set the beam curvature to 0

    for beam_num = 1:num_beams

        %Find displacement of current beam that minimizes substrate energy
        subs_disp = sum(eq_pos(beam_num:end)-beam_ledge_base(beam_num:end))/(num_beams+1-beam_num);

        %Find displacement to left (closer to symmetry) neighboring beam
        if beam_num > 1
            gapLR = beam_ledge_base(beam_num-1)+pitch-beam_ledge_tip(beam_num); 
        else
            gapLR = -gap/2;
        end

        if gapLR > subs_disp 
            %If |displacement| to previous beam is smaller than displacement to eq_pos AND disp to eq_pos is negative

            beam_disp = gapLR; %Beam will move until it contacts next beam (distance of gapLR)

            %Derivative of energy wrt to position is:
            dU_vec = [q_eff/2*(eq_pos(beam_num)^2-2*eq_pos(beam_num)*(4*beam_disp/15)+((beam_disp/3)^2*88/63)),0,0,0,-16*11*EI*beam_disp^2/5];

            L_eff_sol = roots(dU_vec);

            %find the real solution of where the beam contacts the neighboring beam with the physically correct value 
            for sol_num=1:length(L_eff_sol)
                if isreal(L_eff_sol(sol_num)) && L_eff_sol(sol_num) > 0 && L_eff_sol(sol_num) < width
                            L_eff(beam_num) = L_eff_sol(sol_num); %Set effective length to point of contact
                end
            end

            tip_disp = beam_disp; %Link tip displacement to beam displacement

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

    end
end






function serpentine_res = area_to_resistance(num_beams,b,Rs,Rl,Rc,L_eff)
    %num_beams,Rs,Rl,rho_c,xi_int,b,contact_area
    L_eff_mirr = [fliplr(L_eff(2:num_beams)) L_eff(1) L_eff(2:num_beams)]; %Mirror the effective length matrix
    contact_length = b- L_eff_mirr; %find contact length between beams
    num_meanders = num_beams; %Mirror the serpentine
    contact_res = Rc./contact_length;
      
    serpentine_res = res_network(Rs,Rl,contact_res,num_meanders);
end

%Use contact resistance array to find total serpentine resistance
function serpentine_res = res_network(Rs,Rl,Rc_array,num_meanders)
%inputs are Rs, Rl, Rc array, and num_meanders. output is total resistance
%for array of contact resistances, find the total serpentine resistance

    num_nodes = (num_meanders-2)*4+5;

    G = zeros(num_nodes);
    off_diag_l = 2; %counter for off diagonal elements of contact resistance
    off_diag_m = 2; %counter for off diagonal elements of contact resistance
    contact_res_index = 1;
    
    %Generate a conductance matrix for a specified number of nodes
    for k =1:num_nodes %Step through columns
        %Diagonal elements
        if k == 1
            G(k,k) = -(1/(Rl+Rs)+1/Rc_array(k));
            G(k+1,k)= 1/(Rl+Rs);
            G(k+2,k) = 1/Rc_array(k); %first lower off-diagonal
            G(k,k+2) = 1/Rc_array(k); %first upper off_diagonal
        elseif k == 2
            G(k-1,k) = 1/(Rl+Rs);
            G(k,k) = -(1/(Rl+Rs)+1/Rc_array(k)+1/Rl); %%%%%%
            G(k+1,k) = 1/Rl;
        elseif k == num_nodes
            G(k-1,k) =1/Rl;
            G(k,k) = -(1/(Rl+Rs)+1/Rc_array(contact_res_index)+1/Rl); %%%%%
        else
            if mod(k,2) == 1
                G(k,k) = -(1/Rs+1/Rc_array(contact_res_index)+1/Rl); %Contact resistance changes even/odd by 2
            end
            if mod(k,2) == 0
                G(k,k) = -(1/Rs+1/Rc_array(contact_res_index+2)+1/Rl);
                contact_res_index = contact_res_index + 1;
                %Link number of nodes and meander count
            end
            G(k-1,k) = (1/Rl)*mod(k,2)+(1/Rs)*(1-mod(k,2)); %1/Rl and 1/Rs terms alternate based on matrix column
            G(k+1,k) = (1/Rl)*(1-mod(k,2))+(1/Rs)*mod(k,2);
        end
        if k > 4 && mod(k,2) == 1
            G(k-3,k) = 1/Rc_array(off_diag_l);
            off_diag_l = off_diag_l + 1;
            G(k,k-3) = 1/Rc_array(off_diag_m);
            off_diag_m = off_diag_m + 1;
        end
    end

    %Assume an input test current of 1 A
    I_test = zeros(num_nodes,1);
    I_test(1) = -1;

    %Solve for node voltages using conductance matrix and current
    node_voltages = G\I_test;
    serpentine_res = node_voltages(1);
end











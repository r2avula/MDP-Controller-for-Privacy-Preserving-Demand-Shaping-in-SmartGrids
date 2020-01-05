function [policy_all_partitions] = getMDP_deviation_privacy_tradeoff_policy(evalCacheParams)
%% Hypothesis-aware EMU optimal strategy design
controller_Params = evalCacheParams.controller_Params;
adversary_Params = evalCacheParams.adversary_Params;
controller_essParams_all_partitions = evalCacheParams.controller_essParams;
controller_HMM_params = evalCacheParams.controller_HMM_params;
adversary_HMM_params = evalCacheParams.adversary_HMM_params;
deglifePartitions_num = controller_Params.deglifePartitions_num;
batteryRatedCapacityInAh = controller_Params.batteryRatedCapacityInAh;
paramsPrecisionDigits = controller_Params.paramsPrecisionDigits;

tradeOff_omega_withinUtility = evalCacheParams.tradeOff_omega_withinUtility;
tradeOff_sigma_forcost = evalCacheParams.tradeOff_sigma_forcost;
num_tradeoffs = length(tradeOff_omega_withinUtility);
tradeOff_sigma_forcost_idx = num_tradeoffs;

policy_all_partitions = cell(deglifePartitions_num,num_tradeoffs);

policyParams = struct;
policyParams.controller_Params = controller_Params;
policyParams.controller_HMM_params = controller_HMM_params;
policyParams.adversary_Params = adversary_Params;
policyParams.adversary_HMM_params = adversary_HMM_params;
for partition_idx = 1:deglifePartitions_num
    for tradeOff_omega_withinUtility_idx = 2:num_tradeoffs-1
            controller_essParams = controller_essParams_all_partitions{partition_idx};
            policyParams.controller_essParams = controller_essParams;
            
            policyParams.privacyWeight = round(tradeOff_sigma_forcost(tradeOff_sigma_forcost_idx)*(1 - tradeOff_omega_withinUtility(tradeOff_omega_withinUtility_idx))*evalCacheParams.privacyWeight,paramsPrecisionDigits);
            policyParams.deviationWeight = round(tradeOff_sigma_forcost(tradeOff_sigma_forcost_idx)*tradeOff_omega_withinUtility(tradeOff_omega_withinUtility_idx)*evalCacheParams.deviationWeight,paramsPrecisionDigits);
                        
            policyFileNamePrefix = strcat('cache/policy_',...
                num2str(batteryRatedCapacityInAh,'%02d'),'_',...
                num2str(partition_idx),'_',...
               num2str(tradeOff_omega_withinUtility_idx),'_',...
               num2str(tradeOff_sigma_forcost_idx),'_');
            
            [policyFileName,fileExists] = findFileName(policyParams,policyFileNamePrefix,'policyParams');
            if(fileExists)
                disp(strcat({'MDP control policy found in: '},policyFileName));
                load(policyFileName,'policy');
            else
                k_num_in_day = controller_Params.k_num_in_day;
                h_num = controller_Params.h_num;
                controller_x_num = controller_Params.x_num;
                controller_y_num = controller_Params.y_num;
                z_num = controller_Params.z_num;
                d_num = controller_Params.d_num;
                p_pu = controller_Params.p_pu;
                
                timeHorizonsPerDay = controller_Params.timeHorizonsPerDay;
                
                k_num_in_horizon = k_num_in_day/timeHorizonsPerDay;
                
                d_offset = controller_Params.d_offset;
                x_offset = controller_Params.x_offset;
                y_offset = controller_Params.y_offset;
                C_HgHh = controller_Params.C_HgHh;                
                
                privacyWeight = policyParams.privacyWeight;
                deviationWeight = policyParams.deviationWeight;
                
                paramsPrecisionDigits = controller_Params.paramsPrecisionDigits;
                beliefSpacePrecisionDigits = controller_Params.beliefSpacePrecisionDigits;
                paramsPrecision = 10^(-paramsPrecisionDigits);
                
                controller_M_b = controller_HMM_params.M_b;
                
                adversary_P_Hk = adversary_HMM_params.P_Hk;
                adversary_P_HgHn1 = adversary_HMM_params.P_HgHn1;
                adversary_M_b = adversary_HMM_params.M_b;
                adversary_x_num = adversary_Params.x_num;
                
                
                z_kp1_idx_map = controller_essParams.z_kp1_idxs_map;
                
                deltaVectors = eye(h_num);
                
                controllerDecisions = cell(k_num_in_day,1);
                priorBeliefSpacePartitionFirstElements = cell(k_num_in_day+1,1);
                
                priorBeliefSpacePartitionFirstElements{k_num_in_day+1} = unique([0;1]);
                priorBeliefSpacePartition_num_Np1 = length(priorBeliefSpacePartitionFirstElements{k_num_in_day+1})-1;
                valueFuncation_kp1 = zeros(z_num,controller_x_num,h_num,controller_y_num,priorBeliefSpacePartition_num_Np1);
                
                U_delta_hs_num = h_num*(h_num-1);
                U_delta_hs = zeros(U_delta_hs_num,h_num);
                U_delta_vertices_first_elem = zeros(U_delta_hs_num,1);
                I_h_num = eye(h_num);
                idx = 1;
                for i_idx = 1:U_delta_hs_num
                    for j_idx = 1:U_delta_hs_num
                        if(j_idx ~=i_idx)
                            U_delta_hs(idx,:) = (I_h_num(i_idx,:)-I_h_num(j_idx,:))*C_HgHh;
                            idx = idx + 1;
                        end
                    end
                end
                for hp_idx = 1:U_delta_hs_num
                    U_delta_vertices_first_elem(hp_idx) = U_delta_hs(hp_idx,2)/(U_delta_hs(hp_idx,2) - U_delta_hs(hp_idx,1));
                end
                
                U_delta_vertices_first_elem(U_delta_vertices_first_elem<0) = [];
                U_delta_vertices_first_elem = unique([0;round(U_delta_vertices_first_elem,beliefSpacePrecisionDigits);1]);
                
                disp('Computing policy...');
                for horizon_idx = timeHorizonsPerDay:-1:1
                    for k_in_horizon = k_num_in_horizon:-1:1
                        k_in_day = (horizon_idx-1)*k_num_in_horizon + k_in_horizon;
                        
                        J_kp1_vertices_first_elem = unique([0;priorBeliefSpacePartitionFirstElements{k_in_day+1};U_delta_vertices_first_elem;1]);
                        J_kp1_vertices_num = length(J_kp1_vertices_first_elem);
                        
                        J_kp1_inv_vertices_num = J_kp1_vertices_num*(adversary_x_num+1);
                        J_kp1_inv_vertices = zeros(h_num,J_kp1_inv_vertices_num);
                        
                        for y_k_idx = 1:adversary_x_num
                            if(abs(det(adversary_M_b(:,:,y_k_idx,k_in_day)))>=paramsPrecision)
                                temp = adversary_M_b(:,:,y_k_idx,k_in_day)\transpose([J_kp1_vertices_first_elem,1-J_kp1_vertices_first_elem]);
                                J_kp1_inv_vertices(:,(y_k_idx-1)*J_kp1_vertices_num+1:y_k_idx*J_kp1_vertices_num) = temp;
                            end
                        end
                        if(abs(det(adversary_P_HgHn1(:,:,k_in_day)))>=paramsPrecision)
                            temp = adversary_P_HgHn1(:,:,k_in_day)\transpose([J_kp1_vertices_first_elem,1-J_kp1_vertices_first_elem]);
                            J_kp1_inv_vertices(:,(adversary_x_num)*J_kp1_vertices_num+1:(adversary_x_num+1)*J_kp1_vertices_num) = temp;
                        end
                        
                        idx = zeros(J_kp1_inv_vertices_num,1);
                        for h_idx = 1:h_num
                            idx(J_kp1_inv_vertices(h_idx,:)<0)=1;
                        end
                        J_kp1_inv_vertices(:,idx==1) = [];
                        J_kp1_inv_vertices(:,sum(J_kp1_inv_vertices,1)==0)= [];
                        J_kp1_inv_vertices_num = size(J_kp1_inv_vertices,2);
                        for idx = 1:J_kp1_inv_vertices_num
                            J_kp1_inv_vertices(:,idx) = J_kp1_inv_vertices(:,idx)/sum(J_kp1_inv_vertices(:,idx));
                        end
                        
                        G_k_vertices_first_elem_t = transpose(J_kp1_inv_vertices(1,:));
                        
                        fs_filter_vertices_first_elem = zeros(adversary_x_num,1);
                        for y_k_idx = 1:adversary_x_num
                            fs_filter_hp_coeffs = ones(1,h_num)*adversary_M_b(:,:,y_k_idx,k_in_day);
                            fs_filter_vertices_first_elem(y_k_idx) = (paramsPrecision - fs_filter_hp_coeffs(2))/(fs_filter_hp_coeffs(1)-fs_filter_hp_coeffs(2));
                        end
                        
                        beliefSpacePartitionFirstElements_k = [G_k_vertices_first_elem_t;fs_filter_vertices_first_elem];
                        beliefSpacePartitionFirstElements_k(beliefSpacePartitionFirstElements_k<0|beliefSpacePartitionFirstElements_k>1|isnan(beliefSpacePartitionFirstElements_k)|isinf(beliefSpacePartitionFirstElements_k)) = [];
                        beliefSpacePartitionFirstElements_k = unique([0;round(beliefSpacePartitionFirstElements_k,beliefSpacePrecisionDigits);1]);
                        beliefSpacePartition_num_k = length(beliefSpacePartitionFirstElements_k)-1;
                        
                        partitionTransitionMap_k = zeros(adversary_x_num,beliefSpacePartition_num_k);
                        optAdvStratVectorIdx_k = zeros(adversary_x_num,beliefSpacePartition_num_k);
                        for p_idx = 1:beliefSpacePartition_num_k
                            partition_vertices_first_elem = [beliefSpacePartitionFirstElements_k(p_idx);beliefSpacePartitionFirstElements_k(p_idx+1)];
                            centroid_first_elem = mean(partition_vertices_first_elem);
                            pi_hat_kn1 = round([centroid_first_elem;1-centroid_first_elem],paramsPrecisionDigits);
                            for y_k_idx = 1:adversary_x_num
                                pi_hat_k = adversary_M_b(:,:,y_k_idx,k_in_day)*pi_hat_kn1;
                                if(sum(pi_hat_k)>=paramsPrecision)
                                    pi_hat_k = pi_hat_k/sum(pi_hat_k);
                                else
                                    pi_hat_k = adversary_P_HgHn1(:,:,k_in_day)*pi_hat_kn1;
                                end
                                pi_hat_k = round(pi_hat_k,beliefSpacePrecisionDigits);
                                p_idx_kp1 = find(priorBeliefSpacePartitionFirstElements{k_in_day+1}>=pi_hat_k(1),1)-1;
                                if(p_idx_kp1==0)
                                    p_idx_kp1 = 1;
                                end
                                partitionTransitionMap_k(y_k_idx,p_idx) = p_idx_kp1;
                                risk_fn_k = round(deltaVectors*C_HgHh*pi_hat_k,paramsPrecisionDigits);
                                optAdvStratVectorIdx_k(y_k_idx,p_idx) = find(risk_fn_k == max(risk_fn_k),1);
                            end
                        end
                        
                        controllerDecisions_k = zeros(z_num,controller_x_num,h_num,controller_y_num,beliefSpacePartition_num_k);
                        valueFunction_k = zeros(z_num,controller_x_num,h_num,controller_y_num,beliefSpacePartition_num_k);
                        opt_start = tic;
                        
                        for p_idx = 1:beliefSpacePartition_num_k                            
                            for y_kn1_idx = 1:controller_y_num
                                for h_k_idx = 1:h_num
                                    for x_k_idx = 1:controller_x_num
                                        for z_k_idx = 1:z_num
                                            alphaVector_k = zeros(d_num,1);
                                            for d_k_idx = 1:d_num
                                                z_kp1_idx_t = z_kp1_idx_map(z_k_idx,d_k_idx);
                                                if(isnan(z_kp1_idx_t)||z_kp1_idx_t<=0)
                                                    alphaVector_k(d_k_idx) = inf;
                                                else
                                                    y_k_idx = min((x_k_idx+x_offset) + (d_k_idx+d_offset) - y_offset,controller_y_num);
                                                    if(y_k_idx<1)
                                                        alphaVector_k(d_k_idx) = inf;
                                                    else
                                                        p_kp1_idx = partitionTransitionMap_k(min(y_k_idx,adversary_x_num),p_idx);
                                                        H_guess_k = optAdvStratVectorIdx_k(min(y_k_idx,adversary_x_num),p_idx);
                                                        temp = privacyWeight*C_HgHh(h_k_idx,H_guess_k) + deviationWeight*abs((y_kn1_idx - y_k_idx)*p_pu);
                                                        if(k_in_day<k_num_in_day)
                                                            for x_kp1_idx = 1:controller_x_num
                                                                for h_kp1_idx = 1:h_num
                                                                    temp = temp + valueFuncation_kp1(z_kp1_idx_t,x_kp1_idx,h_kp1_idx,y_k_idx,p_kp1_idx)*...
                                                                        controller_M_b(h_kp1_idx,h_k_idx,x_kp1_idx,k_in_day+1);
                                                                end
                                                            end
                                                        end
                                                        alphaVector_k(d_k_idx) = temp;
                                                    end
                                                end
                                            end
                                            
                                            [opt_val,d_k_idx_star] = min(alphaVector_k);
                                            if(d_k_idx_star<1||d_k_idx_star>d_num)
                                                error('Something is wrong! Wrong "d_k_idx_star" obtained during optimization!')
                                            end
                                            controllerDecisions_k(z_k_idx,x_k_idx,h_k_idx,y_kn1_idx,p_idx) = d_k_idx_star;
                                            valueFunction_k(z_k_idx,x_k_idx,h_k_idx,y_kn1_idx,p_idx) = opt_val;
                                        end
                                    end
                                end
                            end
                        end
                        
                        opt_time = toc(opt_start);
                        Postprocess_start = tic;
                        
                        if(h_num==2) % only for 2D case. Combining in higher dimensions can lead to nonconvex paritions.
                            if(beliefSpacePartition_num_k>1)
                                rows_count= 2*controller_x_num*h_num*z_num*controller_y_num;
                                p_sim_check_matrix = zeros(beliefSpacePartition_num_k,rows_count);
                                for p_idx = 1:beliefSpacePartition_num_k
                                    offset = 0;
                                    for y_kn1_idx = 1:controller_y_num
                                        for z_k_idx = 1:z_num
                                            for x_k_idx = 1:controller_x_num
                                                for h_k_idx = 1:h_num
                                                    range = offset + 1: offset+ 2;
                                                    offset = offset + 2;
                                                    p_sim_check_matrix(p_idx,range) = [controllerDecisions_k(z_k_idx,x_k_idx,h_k_idx,y_kn1_idx,p_idx);valueFunction_k(z_k_idx,x_k_idx,h_k_idx,y_kn1_idx,p_idx)];
                                                end
                                            end
                                        end
                                    end
                                end
                                
                                [~,temp,idxs_in_unique_matrix_corresponding_to_originals_vec] = unique(p_sim_check_matrix,'rows');
                                
                                if(length(temp)~=beliefSpacePartition_num_k)
                                    similarity_partitions_idxs = zeros(beliefSpacePartition_num_k,2);
                                    for p_idx = 1:beliefSpacePartition_num_k
                                        B = idxs_in_unique_matrix_corresponding_to_originals_vec(p_idx+1:end);
                                        q_idx = find(B~=idxs_in_unique_matrix_corresponding_to_originals_vec(p_idx),1);
                                        similarity_partitions_idxs(p_idx,1) = p_idx;
                                        if(isempty(q_idx))
                                            similarity_partitions_idxs(p_idx,2) = p_idx;
                                        else
                                            similarity_partitions_idxs(p_idx,2) = p_idx + q_idx-1;
                                        end
                                    end
                                    
                                    similarity_partitions_end_idxs = unique(similarity_partitions_idxs(:,2));
                                    similarity_partitions_start_idxs = similarity_partitions_end_idxs+1;
                                    similarity_partitions_start_idxs = [1;similarity_partitions_start_idxs(1:end-1)];
                                    beliefSpacePartition_with_similarity_num_k = length(similarity_partitions_start_idxs);
                                    beliefSpacePartitionFirstElements_with_similarity_k = ones(beliefSpacePartition_with_similarity_num_k+1,1);
                                    for p_idx = 1:beliefSpacePartition_with_similarity_num_k
                                        beliefSpacePartitionFirstElements_with_similarity_k(p_idx) = beliefSpacePartitionFirstElements_k(similarity_partitions_start_idxs(p_idx));
                                    end
                                    
                                    controllerDecisions_with_similarity_k = zeros(z_num,controller_x_num,h_num,controller_y_num,beliefSpacePartition_with_similarity_num_k);
                                    controllerDecisions_with_similarity_k(:,:,:,:,1:beliefSpacePartition_with_similarity_num_k) = controllerDecisions_k(:,:,:,:,similarity_partitions_start_idxs);
                                    valueFuncation_with_similarity_k = zeros(z_num,controller_x_num,h_num,controller_y_num,beliefSpacePartition_with_similarity_num_k);
                                    valueFuncation_with_similarity_k(:,:,:,:,1:beliefSpacePartition_with_similarity_num_k) = valueFunction_k(:,:,:,:,similarity_partitions_start_idxs);                                    
                                    wrong_indices_count = sum(controllerDecisions_with_similarity_k(:)<1 | controllerDecisions_with_similarity_k(:) > d_num);
                                    if(wrong_indices_count>0)
                                        error('Something is wrong! Some wrong decisions found after processing the control decisions!');
                                    end
                                    controllerDecisions{k_in_day} = controllerDecisions_with_similarity_k;
                                    valueFuncation_kp1 = valueFuncation_with_similarity_k;
                                    priorBeliefSpacePartitionFirstElements{k_in_day} = beliefSpacePartitionFirstElements_with_similarity_k;
                                else
                                    controllerDecisions{k_in_day} = controllerDecisions_k;
                                    valueFuncation_kp1 = valueFunction_k;
                                    priorBeliefSpacePartitionFirstElements{k_in_day} = beliefSpacePartitionFirstElements_k;
                                end
                            else
                                controllerDecisions{k_in_day} = controllerDecisions_k;
                                valueFuncation_kp1 = valueFunction_k;
                                priorBeliefSpacePartitionFirstElements{k_in_day} = beliefSpacePartitionFirstElements_k;
                            end
                        else
                            controllerDecisions{k_in_day} = controllerDecisions_k;
                            valueFuncation_kp1 = valueFunction_k;
                            priorBeliefSpacePartitionFirstElements{k_in_day} = beliefSpacePartitionFirstElements_k;
                        end
                        Postprocess_time = toc(Postprocess_start);
                        disp(strcat('Time index: ',num2str(k_in_day),...
                            ', total partitions: ',num2str(beliefSpacePartition_num_k), ', unique partitions:',num2str(length(priorBeliefSpacePartitionFirstElements{k_in_day})-1),...
                            ', Opt:',num2str(opt_time),', Postproc: ',num2str(Postprocess_time)));
                    end
                    
                    if(horizon_idx>1)
                        % for the next horizon
                        valueFuncation_kp1_original = valueFuncation_kp1;
                        priorBeliefSpacePartition_num_Np1 = 1;
                        valueFuncation_Np1 = zeros(z_num,controller_x_num,h_num,controller_y_num,priorBeliefSpacePartition_num_Np1);
                        priorBelief = adversary_P_Hk(:,k_in_day-1);
                        
                        p_idx_N = find(priorBeliefSpacePartitionFirstElements{k_in_day}>=priorBelief(1),1)-1;
                        if(p_idx_N==0)
                            p_idx_N = 1;
                        end
                        
                        valueFuncation_Np1(:,:,:,:,1) = valueFuncation_kp1_original(:,:,:,:,p_idx_N);
                        priorBeliefSpacePartitionFirstElements{k_in_day} = unique([0;1]);
                        valueFuncation_kp1 = valueFuncation_Np1;
                    end
                end
                
                %% save designed policy
                
                policy = struct;
                policy.priorBeliefSpacePartitionFirstElements = priorBeliefSpacePartitionFirstElements;
                policy.controllerDecisions = controllerDecisions;
                
                save(policyFileName,'policy','policyParams')
                disp(strcat({'MDP control policy saved in: '},policyFileName));
            end
            policy_all_partitions{partition_idx,tradeOff_omega_withinUtility_idx} = policy;
    end    
end
end



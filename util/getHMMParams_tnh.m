function [HMM_params] = getHMMParams_tnh(params,trainingSMdata,trainingGTdata)
x_num = params.x_num;
x_offset = params.x_offset;
p_pu = params.p_pu;
h_num = params.h_num;
k_num_in_day = params.k_num_in_day;
beliefSpacePrecisionDigits = params.beliefSpacePrecisionDigits;

P_Hk = zeros(h_num,k_num_in_day);
P_HgHn1 = zeros(h_num,h_num,k_num_in_day);
P_XgH = zeros(x_num,h_num,k_num_in_day);
M_b = zeros(h_num,h_num,x_num,k_num_in_day);

x_k_idxs = min(round(trainingSMdata/p_pu)-x_offset,x_num);

for k_low_res = 1:k_num_in_day
    for h_idx = 1:h_num
        P_Hk(h_idx,k_low_res) = sum(trainingGTdata(k_low_res,:)==h_idx);
    end
    temp = sum(P_Hk(:,k_low_res));
    if(temp >0)
        P_Hk(:,k_low_res) = P_Hk(:,k_low_res)/temp;
    else
        P_Hk(:,k_low_ress) = 1/h_num;
    end   
    
    if(k_low_res==1)
        prev_h_state = trainingGTdata(k_num_in_day,:);
    else
        prev_h_state = trainingGTdata(k_low_res-1,:);
    end
    cur_h_state = trainingGTdata(k_low_res,:);
    
    for hn1_idx = 1:h_num
        for h_idx = 1:h_num
            P_HgHn1(h_idx,hn1_idx,k_low_res) = sum(prev_h_state==hn1_idx & cur_h_state==h_idx);
        end
        if(sum(P_HgHn1(:,hn1_idx,k_low_res))>0)
            P_HgHn1(:,hn1_idx,k_low_res) = P_HgHn1(:,hn1_idx,k_low_res)/sum(P_HgHn1(:,hn1_idx,k_low_res));
        else
            P_HgHn1(:,hn1_idx,k_low_res) = 1/h_num;
        end
    end    
    
    cur_x_state = x_k_idxs(k_low_res,:);
    
    for h_idx = 1:h_num
        for x_idx = 1:x_num
            P_XgH(x_idx,h_idx,k_low_res) = sum(cur_x_state==x_idx & cur_h_state==h_idx);
        end
        if(sum(P_XgH(:,h_idx,k_low_res)) >0)
            P_XgH(:,h_idx,k_low_res) = P_XgH(:,h_idx,k_low_res)/sum(P_XgH(:,h_idx,k_low_res));
        else
            P_XgH(:,h_idx,k_low_res) = 1/h_num;
        end
    end    
end

P_Hk = round(P_Hk,beliefSpacePrecisionDigits);
P_HgHn1 = round(P_HgHn1,beliefSpacePrecisionDigits);
P_XgH = round(P_XgH,beliefSpacePrecisionDigits);

for k_low_res = 1:k_num_in_day
    for x_k_idx=1:x_num
        for h_kn1_idx=1:h_num
            for h_k_idx = 1:h_num
               M_b(h_k_idx,h_kn1_idx,x_k_idx,k_low_res) = P_HgHn1(h_k_idx,h_kn1_idx,k_low_res)*P_XgH(x_k_idx,h_k_idx,k_low_res);
            end
        end
    end
end

M_b = round(M_b,beliefSpacePrecisionDigits);

HMM_params = struct;
HMM_params.P_Hk = P_Hk;
HMM_params.P_HgHn1 = P_HgHn1;
HMM_params.P_XgH = P_XgH;
HMM_params.M_b = M_b;
end

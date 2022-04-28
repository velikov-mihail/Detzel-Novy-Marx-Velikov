clear
clc

% Load dates vector & factors
load dates
load ret

% Original factors
load ff
load hxz
load hmlm 

% Replications
load ff_rep
load hxz_rep
load hmlm_rep

% Extended factors
load ff_c
load frk

% Mitigated facotrs
load umd_sS
load hxz_sS
load hmlm_sS

% Trading costs
load ff_tc_dnmv
load hxz_tc
load hmlm_tc
load frk_tc
load umd_sS_tc
load hxz_sS_tc
load hmlm_sS_tc

% dWs
load ff_dW 
load hxz_dW
load hmlm_dW
load ff_c_dW
load umd_sS_dW 
load hxz_sS_dW
load hmlm_sS_dW

% Main factors structure
factor_struct_main = struct;
factor_list = {'mkt','smb','hml','umd','rmw','cma','rme','rroe','ria','hmlm','rmw_c'};
factor_struct_main(1).label = 'mkt';
factor_struct_main(1).gross_factor = mkt;
factor_struct_main(1).net_factor = mkt;
factor_struct_main(1).tc = zeros(size(mkt));
factor_struct_main(1).to = zeros(size(mkt));
factor_struct_main(1).dW = zeros(size(ret));

for i=2:length(factor_list)
    factor_struct_main(i).label = factor_list(i);
    eval(['factor_struct_main(i).gross_factor = ',char(factor_list(i)),';']);
    eval(['factor_struct_main(i).net_factor = ',char(factor_list(i)),'-',char(factor_list(i)),'_tc;']);
    factor_struct_main(i).tc = factor_struct_main(i).gross_factor-factor_struct_main(i).net_factor;
    eval(['factor_struct_main(i).to = mean(',char(factor_list(i)),'_TO,2);']);            
    eval(['factor_struct_main(i).dW = ',char(factor_list(i)),'_dW;']);            
end

% Replicate the market factor - Ken French's mkt uses returns not adjusted
% for delisting
load ret_x_dl
load me

laggedMe = lag(me, 1, nan);
indIsFiniteRet = isfinite(ret_x_dl);
sumLaggedMe = sum(laggedMe .* indIsFiniteRet,2,'omitnan');
mkt_rep=sum(laggedMe .* ret_x_dl, 2, 'omitnan') ./ sumLaggedMe - rf;

% Replicated factors structure
factor_struct_rep = struct;
factor_list_rep = {'mkt_rep','smb_rep','hml_rep','umd_rep','rmw_rep','cma_rep','rme_rep','rroe_rep','ria_rep','hmlm_rep'};
factor_struct_rep(1).label = 'mkt_rep';
factor_struct_rep(1).gross_factor = mkt_rep;
factor_struct_rep(1).net_factor = mkt_rep;
factor_struct_rep(1).tc = zeros(size(dates));
factor_struct_rep(1).to = nan(size(dates));
factor_struct_rep(1).dW = zeros(size(mkt));

for i=2:length(factor_list_rep)
    factor_struct_rep(i).label = factor_list_rep(i);
    eval(['factor_struct_rep(i).gross_factor = ',char(factor_list_rep(i)),';']);
    eval(['factor_struct_rep(i).net_factor = ',char(factor_list_rep(i)),'-',char(factor_list(i+1)),'_tc;']);
    factor_struct_rep(i).tc = factor_struct_main(i).tc;
    factor_struct_rep(i).to = factor_struct_main(i).to;
    factor_struct_rep(i).dW = factor_struct_main(i).dW;
end

% Mitigated factors structure
factor_struct_sS = struct;
factor_list_sS = {'umd_sS','rme_sS','rroe_sS','ria_sS','hmlm_sS'};

for i=1:length(factor_list_sS)
    factor_struct_sS(i).label = factor_list_sS(i);
    eval(['factor_struct_sS(i).gross_factor = ',char(factor_list_sS(i)),';']);
    eval(['factor_struct_sS(i).net_factor = ',char(factor_list_sS(i)),'-',char(factor_list_sS(i)),'_tc;']);
    factor_struct_sS(i).tc = factor_struct_sS(i).gross_factor-factor_struct_sS(i).net_factor;
    eval(['factor_struct_sS(i).to = mean(',char(factor_list_sS(i)),'_TO,2);']);            
    eval(['factor_struct_sS(i).dW = ',char(factor_list_sS(i)),'_dW;']);            
end


% Freak factor structure
factor_struct_frk(1).label = 'frk';
factor_struct_frk(1).gross_factor = frk;
factor_struct_frk(1).net_factor = frk-frk_tc;
factor_struct_frk(1).tc = frk_tc;
factor_struct_frk(1).to = frk_TO;
factor_struct_frk(1).dW = zeros(size(ret));

% Combine all factors in one big structure
factor_struct = [factor_struct_main factor_struct_rep factor_struct_sS factor_struct_frk];

% Factor model definitions structure
factor_model_defs = struct;
factor_model_defs(1).label = {'FF5'};
factor_model_defs(1).factors = {'mkt','smb','hml','rmw','cma'};
factor_model_defs(2).label = {'FF6'};
factor_model_defs(2).factors = {'mkt','smb','hml','rmw','cma','umd'};
factor_model_defs(3).label = {'HXZ4'};
factor_model_defs(3).factors = {'mkt','rme','rroe','ria'};
factor_model_defs(4).label = {'BS6'};
factor_model_defs(4).factors = {'mkt','smb','rroe','ria','umd','hmlm'};
factor_model_defs(5).label = {'FF5_c'};
factor_model_defs(5).factors = {'mkt','smb','hml','rmw_c','cma'};
factor_model_defs(6).label = {'FF6_c'};
factor_model_defs(6).factors = {'mkt','smb','hml','rmw_c','cma','umd'};
factor_model_defs(7).label = {'FF5_sS'};
factor_model_defs(7).factors = {'mkt','smb','hml','rmw','cma'};
factor_model_defs(8).label = {'FF6_sS'};
factor_model_defs(8).factors = {'mkt','smb','hml','rmw','cma','umd_sS'};
factor_model_defs(9).label = {'HXZ4_sS'};
factor_model_defs(9).factors = {'mkt','rme_sS','rroe_sS','ria_sS'};
factor_model_defs(10).label = {'BS6_sS'};
factor_model_defs(10).factors = {'mkt','smb','rroe_sS','ria_sS','umd_sS','hmlm_sS'};
factor_model_defs(11).label = {'FF5_c_sS'};
factor_model_defs(11).factors = {'mkt','smb','hml','rmw_c','cma'};
factor_model_defs(12).label = {'FF6_c_sS'};
factor_model_defs(12).factors = {'mkt','smb','hml','rmw_c','cma','umd_sS'};

% Store the factor structures
save Data/factors_dnmv factor_struct factor_model_defs dates -v7.3


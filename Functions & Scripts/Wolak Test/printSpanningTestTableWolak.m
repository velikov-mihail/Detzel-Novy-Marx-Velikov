% function printSpanningTestTableWolak(factor_models,factor_model_defs,factor_struct,startend)

%%
clear
clc

load factors_dnmv

s=find(dates==197201);
e=find(dates==202012);

factor_models={'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'};

const=0.01*ones(size(dates));

model_labels=[factor_model_defs.label]';

gross_res=struct;
net_res=struct;

nfacmod=length(factor_models);


ir2g = nan(nfacmod,nfacmod); %Gross squared "Information ratio"
ir2n = nan(nfacmod,nfacmod); %Net squared "Information ratio"
Pctir2g=nan(nfacmod,nfacmod); %Gross %Delta SR(M1,M0) 
Pctir2n=nan(nfacmod,nfacmod); %Net %Delta SR(M1,M0)

pvalsg=nan(nfacmod,nfacmod); %pvals for gross Spanning tests
pvalsgalt=nan(nfacmod,nfacmod); %pvals for gross Spanning tests
pvalsn=nan(nfacmod,nfacmod); %pvals for net Spanning tests

for i=1:nfacmod
    ri=find(strcmp(model_labels,factor_models(i)));
    factorsi=factor_model_defs(ri).factors;
    indi=find(ismember([factor_struct.label],factorsi));    
    gross_retsi=[factor_struct(indi).gross_factor];
    tcostsi=[factor_struct(indi).tc];
    [gw,gSharpei]=calcMve(gross_retsi(s:e,:));    
    
    [nw,nSharpei]=calcNetMve(gross_retsi(s:e,:)-tcostsi(s:e,:),-gross_retsi(s:e,:)-tcostsi(s:e,:));    
    net_reti=gross_retsi(s:e,:)*nw-tcostsi(s:e,:)*abs(nw); %RHS factors = net factors that get +weight in MVE;
    nw
    
    rhf = (repmat(sign(nw)',size(net_reti))).*gross_retsi(s:e,:)-(repmat(abs(sign(nw))',size(net_reti))).*tcostsi(s:e,:);
    rhfzero=find(sum(abs(rhf))==0);
    rhf(:,rhfzero)=[]; 
    rhf=[const(s:e) rhf];

    for j=1:nfacmod
        if i~=j
            rj=find(strcmp(model_labels,factor_models(j)));
            factorsj=factor_model_defs(rj).factors;

            factoriUj=union(factorsj,factorsi);
            indiUj=find(ismember([factor_struct.label],factoriUj));    
            gross_retsiUj=[factor_struct(indiUj).gross_factor];
            [temp, gSharpeiUj]=calcMve(gross_retsiUj(s:e,:));
            ir2g(i,j)=gSharpeiUj^2-gSharpei^2;
            Pctir2g(i,j)=100*ir2g(i,j)/gSharpei^2;
            
            nF=length(factoriUj);
            opts = optimoptions('fmincon','Display','off');
            GrossRets=gross_retsiUj(s:e,:);
            fun = @(w) -sqrt(12)*mean(GrossRets*w)/std(GrossRets*w);
            try
                [wStar,SRStar]=fmincon(fun,temp,[],[],ones(1,nF),1,zeros(nF,1),[],[],opts);
                ir2galt(i,j)=SRStar^2-gSharpei^2;
                Pctir2galt(i,j)=100*ir2galt(i,j)/gSharpei^2;
            end
            
            tcostsiUj=[factor_struct(indiUj).tc];        
            [nw,nSharpeiUj]=calcNetMve(gross_retsiUj(s:e,:)-tcostsiUj(s:e,:),-gross_retsiUj(s:e,:)-tcostsiUj(s:e,:));
            ir2n(i,j)=nSharpeiUj^2-nSharpei^2;
            Pctir2n(i,j)=100*ir2n(i,j)/nSharpei^2;

            %Spanning
            factorij=setdiff(factorsj,factorsi);
            if size(factorij,2)>0
                i
                j
                %Gross 
                indij=find(ismember([factor_struct.label],factorij));    
                gross_retsij=[factor_struct(indij).gross_factor];
                lhfg = [gross_retsij(s:e,:)];
                rhfg = [const(s:e) gross_retsi(s:e,:)];
                [alphahatg,alphacovg,bg,covbg]=gmmAlphas(lhfg,rhfg);
                pvalsg(i,j)=1-chi2cdf(alphahatg'*inv(alphacovg)*alphahatg,size(alphahatg,1));
                pvalsgalt(i,j)=WolakAlphaTest(alphahatg,alphacovg);
             
                %Net
                tcostsij=[factor_struct(indij).tc];        
                lhf=[gross_retsij(s:e,:)-tcostsij(s:e,:), -gross_retsij(s:e,:)-tcostsij(s:e,:)];
                [alphahat,alphacov,b,covb]=gmmAlphas(lhf,rhf);
                pvalsn(i,j)=WolakAlphaTest(alphahat,alphacov);
            else
            pvalsg(i,j)=1;
            pvalsn(i,j)=1;
            ir2g(i,j)=0;
            Pctir2g(i,j)=0;
            ir2n(i,j)=0;
            Pctir2n(i,j)=0;
            end
        else
            pvalsg(i,j)=nan;
            pvalsn(i,j)=nan;
            ir2g(i,j)=nan;
            Pctir2g(i,j)=nan;
            ir2n(i,j)=nan;
            Pctir2n(i,j)=nan;
        end
    end
end

pvalsn

h=factor_models';
h=regexprep(h,'_c','\\textsubscript{C}');
h=regexprep(h,'_sS1','');

mat2Tex([ir2g nan(size(ir2g,1),1) ir2n],[pvalsg nan(size(pvalsg,1),1) pvalsn],h,3);
mat2Tex([Pctir2g nan(size(ir2g,1),1) Pctir2n],[Pctir2g nan(size(ir2g,1),1) Pctir2n],h,1);

mat2Tex([ir2galt nan(size(ir2galt,1),1) ir2n],[pvalsgalt nan(size(pvalsgalt,1),1) pvalsn],h,3);
mat2Tex([Pctir2galt nan(size(ir2g,1),1) Pctir2n],[Pctir2galt nan(size(ir2g,1),1) Pctir2n],h,1);

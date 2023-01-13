function [Final_Ko_vs_Ct, Final_Ko_vs_Ct_Thr, MetName, FC_Ko_vs_Ct] = EndoIFT88MEF_KoVsCt(Cond)
%
% Metabolomic Analysis of Ko vs Ctrl used in paper "Primary Cilia Respond to Glutamine Availability via Asparagine Synthetase", 
% M. E Steidl, E. A. Nigro, et al., Nature Metabolism, 2023.
%
% Input:
%
% a) Cond = string that allows to compute the results. If Cond = 'HBSS24hs' the function analyzes metabolites in 
% KO vs CT cultured in HBSS for 24 hours. If Cond = 'HBSSQ24hs' the function analyzes metabolites in 
% KO vs CT cultured in HBSS + L-Glutamine for 24 hours.  
%
% Output:
%
% a) Final_Ko_vs_Ct = it containes the list of metabolites having a
% statistical significant fold-change (FDR <= 0.05. FDR values have been derived from Benjamini-Hochberg FDR correction procedure 
% applied to p-values computed by means of Student's t-test). First column
% containes the metabolite name, second column the log2 fold-change, the third
% column the p-value, and the fourth column the adjusted p-value.
% b) Final_Ko_vs_Ct_Thr = it containes the list of metabolites having a
% statistical significant fold-change (FDR <= 0.05) for which absolute
% value is >= 1.5. First column containes the metabolite name, second column the log2 fold-change, the third
% column the p-value, and the fourth column the adjusted p-value.
% c) MetName = name of all the metabolites.
% d) FC_Ko_vs_Ct_Thr = it containes the fold-change anf the p-value for all
% the metabolites in MetName. 

[num,txt,raw] = xlsread('SourceData.xlsx',2);

[n_txt, m_txt] = size(txt);

MetName = {};

Ct = [];

Ko = [];

for i = 3:m_txt
    
    MetName{i-2,1} = txt{1,i};
    
    %%%%%%%%HBSS24hs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(Cond,'HBSS24hs')
    
        Ct(1:5,i-2) = num(1:5,i-2);
    
        Ko(1:5,i-2) = num(11:15,i-2);
        
    end
    
    %%%%%%%%HBSSQ24hs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(Cond,'HBSSQ24hs')
    
        Ct(1:5,i-2) = num(6:10,i-2);
    
        Ko(1:5,i-2) = num(16:20,i-2);
        
    end
        
end

FC_Ko_vs_Ct = [];

Final_Ko_vs_Ct = {};

for i = 1:m_txt-2

    FC_Ko_vs_Ct(i,1) = mean(Ko(:,i))/mean(Ct(:,i));
    
    [h,p] = ttest2(Ko(:,i),Ct(:,i));
    
    FC_Ko_vs_Ct(i,2) = p;
    
end

FDR = mafdr(FC_Ko_vs_Ct(:,2),'BHFDR','true');

FC_Ko_vs_Ct(:,3) = FDR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Final_Ko_vs_Ct{1,1} = 'Metabolite';

Final_Ko_vs_Ct{1,2} = 'FC';

Final_Ko_vs_Ct{1,3} = 'Log2(FC)';

Final_Ko_vs_Ct{1,4} = 'Adj-pval';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Final_Ko_vs_Ct_Thr{1,1} = 'Metabolite';

Final_Ko_vs_Ct_Thr{1,2} = 'FC';

Final_Ko_vs_Ct_Thr{1,3} = 'Log2(FC)';

Final_Ko_vs_Ct_Thr{1,4} = 'Adj-pval';


for i = 1:m_txt-2
    
    if FDR(i) <= 0.05
        
        Final_Ko_vs_Ct{end+1,1} = MetName{i,1};
        
        Final_Ko_vs_Ct{end,2} = FC_Ko_vs_Ct(i,1);
        
        Final_Ko_vs_Ct{end,3} = log2(FC_Ko_vs_Ct(i,1));
        
        Final_Ko_vs_Ct{end,4} = FDR(i);
        
        if abs(log2(FC_Ko_vs_Ct(i,1))) >= abs(log2(1.5)) 
            
            Final_Ko_vs_Ct_Thr{end+1,1} = MetName{i,1};
        
            Final_Ko_vs_Ct_Thr{end,2} = FC_Ko_vs_Ct(i,1);
        
            Final_Ko_vs_Ct_Thr{end,3} = log2(FC_Ko_vs_Ct(i,1));
        
            Final_Ko_vs_Ct_Thr{end,4} = FDR(i);
            
        end
        
    end
    
end

%%%%Volcano Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

hold on

for i = 1:length(FDR)
    
    if  -log10(FDR(i)) <= -log10(0.05) 
        
         plot(log2(FC_Ko_vs_Ct(i,1)),-log10(FDR(i)), 'o', 'MarkerEdgeColor',[224/255,224/255,224/255], 'MarkerFaceColor',[224/255,224/255,224/255],'MarkerSize',10);
         
    else
        
        if log2(FC_Ko_vs_Ct(i,1)) > 0 
            
            plot(log2(FC_Ko_vs_Ct(i,1)),-log10(FDR(i)), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',10);
            
        end 
        
        if log2(FC_Ko_vs_Ct(i,1)) < 0 
            
            plot(log2(FC_Ko_vs_Ct(i,1)),-log10(FDR(i)), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g','MarkerSize',10);
            
        end
        
    end
     
end

xlabel('Log_{2} Fold Change','fontsize',40,'fontweight','b','color','k');

ylabel('-Log_{10} Adjusted p-value','fontsize',40,'fontweight','b','color','k');

title(strcat('Ko vs Ct-',Cond))

set(gca,'fontsize',30,'LineWidth',2,'fontWeight','bold');

hold off

clustergram([Ct',Ko'],'Standardize','Row','ColumnLabels',{'CT','CT','CT','CT','CT','KO','KO','KO','KO','KO'},'RowLabels',MetName,'Colormap',redbluecmap);
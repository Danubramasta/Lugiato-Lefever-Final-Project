function stability(P,Norm,lamb,no_fig,ls)
id = 0;
for idx=1:length(lamb)
    if abs(lamb(idx))<1e-3
        lamb(idx) = 0;
    end
end
temp1 = [];
temp2 = [];
for idx=1:length(lamb)-1
    if (lamb(idx)>0 && lamb(idx+1)<=0)
        temp1 = [temp1 idx];
    elseif (lamb(idx)<=0 && lamb(idx+1)>0)
        temp2 = [temp2 idx];
    end
end
temp = sort([temp1 temp2]);

if length(temp)>1
    Norm_g{1} = Norm(1:temp(1));
    P_g{1} = P(1:temp(1));
    for idx=2:length(temp)
        Norm_g{idx} = Norm(id+temp(idx-1):temp(idx));
        P_g{idx} = P(id+temp(idx-1):temp(idx));
    end
    Norm_g{idx+1} = Norm(id+temp(end):end);
    P_g{idx+1} = P(id+temp(end):end);
elseif length(temp) == 1
    Norm_g{1} = Norm(1:temp(1));
    P_g{1} = P(1:temp(1));
    Norm_g{2} = Norm(id+temp(1):end);
    P_g{2} = P(id+temp(1):end);
else
    Norm_g = Norm;
    P_g = P;
end

figure(no_fig)
hold on
if length(temp)>=1
    for idx = 1:length(Norm_g)
        if lamb(1)<=0
            if mod(idx,2)==1
                plot((P_g{idx}),Norm_g{idx},'b-','LineWidth',ls)
            else
                plot((P_g{idx}),Norm_g{idx},'r-','LineWidth',ls-1)
            end
        else
            if mod(idx,2)==1
                plot((P_g{idx}),Norm_g{idx},'r-','LineWidth',ls-1)
            else
                plot((P_g{idx}),Norm_g{idx},'b-','LineWidth',ls)
            end
        end
    end
else
    if lamb(1)<=0
        plot((P_g),Norm_g,'b-','LineWidth',ls)
    else
        plot((P_g),Norm_g,'r-','LineWidth',ls-1)
    end
end
% if length(temp)>=1
%     for idx = 1:length(Norm_g)
%         if lamb(1)<=0
%             if mod(idx,2)==1
%                 plot(abs(P_g{idx}),Norm_g{idx},'b.','LineWidth',ls)
%             else
%                 plot(abs(P_g{idx}),Norm_g{idx},'r.','LineWidth',ls)
%             end
%         else
%             if mod(idx,2)==1
%                 plot(abs(P_g{idx}),Norm_g{idx},'r.','LineWidth',ls)
%             else
%                 plot(abs(P_g{idx}),Norm_g{idx},'b.','LineWidth',ls)
%             end
%         end
%     end
% else
%     if lamb(1)<=0
%         plot(abs(P_g),Norm_g,'b.','LineWidth',ls)
%     else
%         plot(abs(P_g),Norm_g,'r.','LineWidth',ls)
%     end
% end
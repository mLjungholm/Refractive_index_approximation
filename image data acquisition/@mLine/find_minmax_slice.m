function find_minmax_slice(this)
edgevals = this.points(this.leftEdgeIndex,:);
mms = 2;
% gas = 10;
sgs = 15;
mm = smoothdata(edgevals,'movmedian',mms);
% ga = smoothdata(mm,'gaussian',gas);
sg = smooth(mm,sgs,'sgolay');

% ga_peak_vals = [min(ga);max(ga)];
sg_peak_vals = [min(sg);max(sg)];
% ga_peak_pos = [find(ga == ga_peak_vals(1)); find(ga == ga_peak_vals(2))];
sg_peak_pos = [find(sg == sg_peak_vals(1)); find(sg == sg_peak_vals(2))];

this.leftPhaseMin = sg_peak_pos(1);
this.leftPhaseMax = sg_peak_pos(2);



% x = (1:1:this.imNums)';
% figure(1)
% hold on
% grid minor

% plot(x,mm,'color',[0.9290 0.6940 0.1250])
% plot(x,ga,'color',[0 0.4470 0.7410].*0.8)
% plot(x,sg,'color',[0.6350 0.0780 0.1840])
% plot(x',edgevals,'k-')
% plot(ga_peak_pos,ga_peak_vals,'ko')
% plot(sg_peak_pos,sg_peak_vals,'ko')
% 
% 
% mmtext = strcat('moving median ',' ',num2str(mms));
% gatext = strcat('Gaussian ',' ',num2str(gas));
% sgtext = strcat('Sgolay ',' ',num2str(sgs));
% legend(mmtext,gatext,sgtext,'raw')

end
function [status] = Truss_thickness_plot2D_colormap(xstar,ystar,ID,A)

nele = size(ID,1);
Amin=min(A);
Amax=max(A);
DeltaA=Amax-Amin;

njet = 12;
colorjet = jet(njet);
%colorjet = turbo(njet);
%colorjet = copper(njet);
scol = zeros(nele,3);
for s=1:nele
   scol(s,:) = interp1( (0:(njet-1))/(njet-1), colorjet , (A(s)-Amin)/DeltaA );
end


figure
colormap(jet(njet))
%colormap(turbo(njet))
%colormap(copper(njet))
colorbar('Ticks',[0:1.0/njet:1],'TickLabels',num2str([Amin:DeltaA/njet:Amax]'))
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [2.2 1.6]);
axis equal
hold on
    for i=1:nele
            plot([xstar(ID(i,1)),xstar(ID(i,2))],[ystar(ID(i,1)),ystar(ID(i,2))],'-','Color',scol(i,:),'LineWidth',2.0)
    end
hold off

status = 0;

end


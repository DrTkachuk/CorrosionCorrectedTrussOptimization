function [status] = Truss_thickness_plot2D(xstar,ystar,ID,A)

nele = size(ID,1);
Amin=min(A);
Amax=max(A);
DeltaA=Amax-Amin;

    njet = 24;
    colorjet = jet(njet);

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [2.2 1.6]);
axis equal
hold on
    for i=1:nele
            plot([xstar(ID(i,1)),xstar(ID(i,2))],[ystar(ID(i,1)),ystar(ID(i,2))],'Color',[0,0.7,(A(i)-Amin)/DeltaA],'LineWidth',2.0)
    end
hold off
%colorbar
status = 0;

end


%%% use this code to plot eigenvectors of the Amundsen sea plot 
%%% code segments from Nurbek Tazhimbetov's 2022 PhD thesis
%%% changes by Nestor Walters, advised by Eric Dunham

%%% given a .mat file containing eigenvalues and eigenvectors of the 
% problem Adq/dt = Bq + forcing from Nurbek's code, plots the eigenvectors
% as 2D+color plots, and lists freq and decay in the title
% NOTE:: be sure to adjust for *i/2pi, depending on whether using raw
% eigenvalues or adjusted for frequency

%%% variables
% W := 2N x K matrix containing K eigenvectors of Adq/dt = Bq + forcing
% L := K eigenvalues of the system stored as vector or diag matrix
% here, 2N is the length of the full [w;phi] vector in Nurbek's code
% and K is the number of eigenvals/vecs that have been solved for
% w0 comes from illapel_amundsen_ice.m and gives the length of the
% displacement-only portion of the vector
% see am0_get_eigs.m



%%% get mesh variables
run illapel_amundsen_ice.m
Size = length(w0);

%%% load desired evec matrix
load am0_LambdaW2knorm_smabs_bvar_tvar.mat

X = W2Knorm(1:Size,:); % these may need to be adjusted
L = Lambda2K;          % to what you call the eval, evec matrices

% now to plot
Suffix = ["-ImRe","-Abs"];
% set s = 1 for imag and real plot, s = 2 for abs plot
s = 2;
for k = 1:5
    % normalize w
    w = X(:,k);
    w = w / max(abs(w));
    
    % set title variables
    ptitle = ['AmunBvarTcon-Smabs',num2str(k)];
    suffix = Suffix(s);
    retitle = [ptitle,'-Real'];
    imtitle = [ptitle,'-Imag'];
    abtitle = [ptitle,'-Abs'];
    stitle = ['Real (Decay): ',num2str(real(L(k))),...
        '  Im (Freq): ',num2str(imag(L(k)))];
    
  
  switch s
      case 1

      %%%%% plot the real and imaginary of w %%%%%%
        figure(k)
        ax(1) = subplot(2,2,1);
        hold on
        cmap
        caxis([-1 1])
        for i = 1:length(g_ice.boundaryGroups.FR)
            free_coord = g.getBoundary(g_ice.boundaryGroups.FR(i));
            p1 = plot3(free_coord(:, 1), free_coord(:, 2), 1000 + 0*free_coord(:, 1), 'black-', 'LineWidth', 1.5);
            uistack(p1, 'top')                        
            hold on
        end
        for i = 1:length(g.boundaryGroups.R)
            reflectN = g.getBoundary(g.boundaryGroups.R(i));
            p2 = plot3(reflectN(:, 1), reflectN(:, 2), 1000 + 0*reflectN(:, 1), 'black-', 'LineWidth', 0.5);
            uistack(p2, 'top')
            hold on
        end
        
        xlabel('$x$ (km)', 'Interpreter', 'latex')
        ylabel('$y$ (km)', 'Interpreter', 'latex')
        xlim([plot_details.left, plot_details.right])
        ylim([plot_details.bottom, plot_details.top])
        
        ph = multiblock.Surface(g, real(w));
        ph.set('LineStyle', 'none');
        ph.set('FaceColor', 'interp');
        ph.CData = real(w);
        ph.ZData = real(w);
        
        c = colorbar;
        c.Label.String = '$w$ in meters';
        c.Label.Interpreter = 'latex';
        c.Label.FontSize = 18;
        
        title(retitle)
        subtitle(stitle)
        
        ax(2) = subplot(2,2,2);
        hold on
        cmap
        caxis([-1 1])
        for i = 1:length(g_ice.boundaryGroups.FR)
            free_coord = g.getBoundary(g_ice.boundaryGroups.FR(i));
            p1 = plot3(free_coord(:, 1), free_coord(:, 2), 1000 + 0*free_coord(:, 1), 'black-', 'LineWidth', 1.5);
            uistack(p1, 'top')                        
            hold on
        end
        for i = 1:length(g.boundaryGroups.R)
            reflectN = g.getBoundary(g.boundaryGroups.R(i));
            p2 = plot3(reflectN(:, 1), reflectN(:, 2), 1000 + 0*reflectN(:, 1), 'black-', 'LineWidth', 0.5);
            uistack(p2, 'top')
            hold on
        end
        
        xlabel('$x$ (km)', 'Interpreter', 'latex')
        ylabel('$y$ (km)', 'Interpreter', 'latex')
        xlim([plot_details.left, plot_details.right])
        ylim([plot_details.bottom, plot_details.top])
        
        ph = multiblock.Surface(g, imag(w));
        ph.set('LineStyle', 'none');
        ph.set('FaceColor', 'interp');
        ph.CData = imag(w);
        ph.ZData = imag(w);
        
        c = colorbar;
        c.Label.String = '$w$ in meters';
        c.Label.Interpreter = 'latex';
        c.Label.FontSize = 18;
        
        title(imtitle)
        subtitle(stitle)
  

%%%% plot the absolute of w %%%%%%
      case 2
        figure(k)
        ax(3) = subplot(1,1,1);
        hold on
        colormap(ax(3),"parula")
        clim([0 1])
        for i = 1:length(g_ice.boundaryGroups.FR)
            free_coord = g.getBoundary(g_ice.boundaryGroups.FR(i));
            p1 = plot3(free_coord(:, 1), free_coord(:, 2), 1000 + 0*free_coord(:, 1), 'black-', 'LineWidth', 1.5);
            uistack(p1, 'top')                        
            hold on
        end
        for i = 1:length(g.boundaryGroups.R)
            reflectN = g.getBoundary(g.boundaryGroups.R(i));
            p2 = plot3(reflectN(:, 1), reflectN(:, 2), 1000 + 0*reflectN(:, 1), 'black-', 'LineWidth', 0.5);
            uistack(p2, 'top')
            hold on
        end
        
        xlabel('$x$ (km)', 'Interpreter', 'latex')
        ylabel('$y$ (km)', 'Interpreter', 'latex')
        xlim([plot_details.left, plot_details.right])
        ylim([plot_details.bottom, plot_details.top])
        
        ph = multiblock.Surface(g, abs(w));
        ph.set('LineStyle', 'none');
        ph.set('FaceColor', 'interp');
        ph.CData = abs(w);
        ph.ZData = abs(w);
        
        c = colorbar;
        c.Label.String = '$w$ in meters';
        c.Label.Interpreter = 'latex';
        c.Label.FontSize = 18;
        
        title(abtitle)
        subtitle(stitle)
  end

    %%%% export whatever you have plotted
    ax = gcf;
    picfilename = [ptitle,char(suffix),'.png'];
    exportgraphics(ax,picfilename);
    
    close all

end

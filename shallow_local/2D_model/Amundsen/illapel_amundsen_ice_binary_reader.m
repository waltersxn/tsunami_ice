%%% given a binary file recorded during the full simulation of
% illapel_amundsen_ice.m (IAMI) 
% this code plays back the binary file as video or stores images files at
% desired intervals
%%% 

%%% to record a binary file, set save_to_binary to 'yes' in IAMI
% note that the binary file name must be changed manually 
% when ice property variables are changed (e.g. BendS) 

%%% original code by Nurbek Tazhimbetov 2022
% changes by Nestor Walters advised by Eric Dunham



run illapel_amundsen_ice.m % need this to recover necessary variables

fileID = fopen('illapel_amundsen_ice_q17_binary.bin', 'r');
Size = n_water;

A = fread(fileID, [2*n_water iter], 'single');
disp(iter)

%%% in our experiments, we tested BendS = 13, 14, 15, 16, 17 
% and replayed the stored binary videos
%%% the files were saved as illapel_amundsen_ice_qxx_binary.bin
% where xx = 13,...,17
%%% some sections are specific to video replay, others to png storage
%%% NOTE: must reset time_in_sec to 43200 in IAMI


%%% adjust:: set q = q of binary bin file being read
% these are used to label pictures and videos 
q = 17;

%%% for pics uncomment these 2 lines
titlename = ['Illapel Amundsen Ice q',num2str(q)];
picfilelabel = ['B_', num2str(q)];

%%% for video uncomment these 4 lines
 %vidfilename = ['Illapel_amundsen_ice',num2str(q)]
 %myVideo = VideoWriter(vidfilename, 'MPEG-4');
 %myVideo.FrameRate = 15;
 %open(myVideo)

figure()
hold on
cmap

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

set(gca,'FontSize',15)
xlabel('$x$ (km)', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('$y$ (km)', 'Interpreter', 'latex', 'FontSize', 18)
xlim([plot_details.left, plot_details.right])
ylim([plot_details.bottom, plot_details.top])

% ith_line = fgetl(fileID);
ph = multiblock.Surface(g, A(1:Size, 1));
% caxis([plot_details.caxis_low, plot_details.caxis_up])
caxis([-0.015, 0.015])

c = colorbar;
c.Label.String = 'Wave height, $w$ (m)';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 20;
ph.set('LineStyle', 'none');
ph.set('FaceColor', 'interp');

hText = text(100, -210, append('$t=$ ', num2str(0), 's'), 'Interpreter', 'latex', 'FontSize', 15);
for i = [1:250:iter, iter]
    disp(i)
        delete(hText)
    hText = text(100, -210, append('$t=$ ', num2str(i*16), 's'), 'Interpreter', 'latex', 'FontSize', 15);
%     ith_line = fgetl(fileID);
    ph.CData = A(1:Size, i);
    ph.ZData = A(1:Size, i);
    drawnow
     %%% for video uncomment these 2 lines
     % frame = getframe(gcf); 
     % writeVideo(myVideo, frame);
    legend([p1, p2], {'Ice shelf front', 'Coastline'}, 'Interpreter', 'latex', 'Location', 'southeast')
    title(titlename)
    %%% for pics uncomment this if statement
    if (i*16) > 8000 % determines the interval of images
        ax = gcf;
        picfilename = [picfilelabel,'_time',num2str(i*16),'.png'];
        exportgraphics(ax, picfilename);
    end
end

%%% uncomment for video
%close(myVideo)

% fclose(fileID);
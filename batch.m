% 批量处理硅片读写
%
for i=1:12
    filename = strcat([num2str(i),'.jpg']);
    img = imread(filename);
    img = img(:,:,1);
    img = imnoise(img,'salt & pepper');
    
    [si,img_out] = Si_counter(img);
    si_loc = (si.end(:)+si.start(:))/2;
    
    [m,n] = size(img_out);
    imshow(img);hold('on');
    for k = 1:length(si.ind)
        % 绘制第k条直线
        plot([si_loc(k),si_loc(k)],[ceil(m/3),floor(2*m/3)],'LineWidth',2,'Color','green');
    end
    text(ceil(n/10),ceil(m/10),strcat(['total num:',num2str(length(si.ind))]),...
        'Color','green','FontSize',24,"FontWeight","bold");
    hold('off');
    ac = gca;
    exportgraphics(ac,strcat(['count_',filename]),'resolution',500,'backgroundcolor',[0,0,0])
    disp(strcat([filename,' saved']))
end
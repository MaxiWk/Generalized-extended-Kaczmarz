% saves all open figures as png files into folder with FolderName

function save_figures(FolderName)

    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = num2str(iFig);
      saveas(FigHandle, fullfile(FolderName, FigName),'png');
    end

end

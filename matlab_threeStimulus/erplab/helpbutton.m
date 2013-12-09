% Author: Javier Lopez-Calderon

try
      [img,map]    = imread('erplab_help.jpg');
      colormap(map)
      [row,column] = size(img);
      p = get(handles.pushbutton_help,'Position');
      w = p(3); h = p(4);
      x = ceil(row/(w*5));
      y = ceil(column/(h*39));
      imgbutton    = img(1:x:end,1:y:end,:);
      set(handles.pushbutton_help,'CData',imgbutton);
catch
      set(handles.pushbutton_help,'String','Help');
end
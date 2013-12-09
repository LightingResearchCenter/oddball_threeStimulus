% Author: Javier Lopez-Calderon

try
      [img,map]    = imread('PDFlogo.jpg');
      colormap(map)
      [row,column] = size(img);
      p = get(handles.pushbutton_pdf,'Position');
      w = p(3); h = p(4);
      x = ceil(row/(w*4.5));
      y = ceil(column/(h*37));
      imgbutton    = img(1:x:end,1:y:end,:);
      set(handles.pushbutton_pdf,'CData',imgbutton);
catch
      set(handles.pushbutton_pdf,'String','Help');
end
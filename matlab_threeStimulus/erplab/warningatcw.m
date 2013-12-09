function warningatcw(message, color)

linewrngn = ['\n' repmat('*',1,100) '\n'];
msgwrng = [linewrngn message linewrngn];
try cprintf(color, '%s', sprintf(msgwrng));catch,fprintf('%s', sprintf(msgwrng));end ;
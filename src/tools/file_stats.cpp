#include "file_stats.h"
#include <stdio.h>
#include "sys/stat.h"
#include "stdlib.h"

int file_stats(char* filename)
{
	struct stat sb;

	if (stat(filename, &sb) == -1)
	{
		return -1;
	}
	
	switch (sb.st_mode & S_IFMT) 
	{
		case S_IFBLK:	return 0;	break;//block device
		case S_IFCHR:	return 0;	break;//character device
		case S_IFDIR:	return 1;	break;//directory
		case S_IFIFO:	return 0;	break;//FIFO/pipe
		case S_IFLNK:	return 2;	break;//symlink
		case S_IFREG:	return 3;	break;//regular file
		case S_IFSOCK:	return 0;	break;//socket
		default:		return 0;	break;
	}
}

bool fexist(char* filename)
{
	return file_stats(filename)>1; // symlink or regular file
}


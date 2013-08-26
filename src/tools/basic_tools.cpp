
#include "tools/basic_tools.h"

std::vector<char*> separate(char* str, char sep)
{
	std::vector<char*> ret;
	int i=0;
	ret.push_back(str);
	while (str[i]!=0)
	{
		if (str[i]==sep)
		{
			str[i]=0;
			ret.push_back(str+i+1);
		}
		i++;
	}
	return ret;
}

std::vector<char*> get_fields(char* line)
{
	std::vector<char*> ret;
	int i=0;
	ret.push_back(line);
	while (line[i]!=0)
	{
		if (line[i]=='\t')
		{
			line[i]=0;
			ret.push_back(line+i+1);
		}
		i++;
	}
	return ret;
}


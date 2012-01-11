/*
This program is parses file formats used for testing Biochemical models
The original file was provided 
Gennemark and Wedelin, Benchmarks for identification of ordinary 
differential equations from time series data,
 Bioinformatics 25(6):780-6

Modified by
Author: Subhodeep Moitra (smoitra@cs.cmu.edu)

BSD License

*/

#include<Python.h>
#include<stdio.h>
#include<stddef.h>
#include<string.h>
#include<stdlib.h>



/*
int main(int argc, char* argv[]) {
	char *inputf = "allProblems/ss_5genes1";
	printf("%s\n",inputf);
	readFile(inputf);
	return 0;
}

*/

PyObject* readFile(const char *fileame);


static PyObject* 
cparser_parse(PyObject* self, PyObject* args){
	const char *fname;
	PyObject *parsedict;
	
	if(!PyArg_ParseTuple(args,"s",&fname))
		return NULL;
	
	parsedict = readFile(fname);
	return parsedict;
}

static PyMethodDef CParserMethods[] = {
	{"parse",cparser_parse, METH_VARARGS,
	"Parse the experiment file and pass it to Python"},
	{NULL,NULL,0,NULL}
};

PyMODINIT_FUNC
initcparser(void) {

	(void) Py_InitModule("cparser",CParserMethods);
}

int
endOfLine(FILE *infile) // checks if we have eol in front of us
{
	char c;
	while((c=getc(infile))==' ');  // eat blanks
	ungetc(c,infile); // put back non-blank character
	if(c=='\n') return 1;
	else return 0;
}


PyObject*
readFile(const char *filename) /* this routine does not change */
{
	FILE *infile;
	char rowbuffer[1000];   /* not so good really */
	char  *keyword, *indexString;
	char keyToken[100];
	int keyIndex,end=0;

	// Init Python objects
	PyObject *parsedict = NULL,*pydict,*pylist;
	PyObject *pykey,*pyval,*pystr,*pyfloat,*pyint;
	PyObject *pyvarlist;
	char *key,*val;
	parsedict = PyDict_New();	
	pyvarlist = PyList_New(0);

	char formatVersion[100], type[100],date[100],url[100],rName[100],lvName[100],name[100];
	char lvName2[100],lpName[100];
	char property[100],reactionType[100],varSetX1[100],rangeK1[100],varSetX2[100],object[100];
	char errorType[100],equation[100];
	int varIndex,index,lvIndex,lvIndex2,lpIndex;
	double lb,ub,lambda,time,value;
	
	char problemName[40];
	char buf[100];
	infile=fopen(filename,"r");
	if(!infile) {printf("Error: file ""%s"" does not exist!\n",filename); return NULL;}

	/* Ununused variables - maybe used later */
	// int lines;
	//int rowlength;
	//char c;
	//int l,w;
	//int accept=0,i;
	//int setIndex ;

	while(1) {
		
		// find the first identifier = keyToken
		
		keyToken[0]='\0';
		if(fscanf(infile,"%s",keyToken)==EOF) break; // find keyToken
		
		// split keyToken into keyword_indexString
		
		keyword=keyToken;
		indexString=keyToken;
		while(*indexString != '\0') {
			if(*indexString == '_') {*indexString='\0'; indexString++; break;}
			indexString++;
		}
		sscanf(indexString,"%d",&keyIndex); // translate indexString to integer

		// now check the keyword and read accordingly
		// in this code we choose to reprint the problem to stdout, this can be replaced as appropriate
			
		if(!strcmp(keyword,"begin")) {
			fscanf(infile," problem %s",problemName);
			printf("begin problem %s\n",problemName);
			pystr= PyString_FromString(problemName);
			PyDict_SetItemString(parsedict,"problemname",pystr);
			end=0;
		}
		else if(!strcmp(keyword,"end")) {
			fscanf(infile," problem %s",problemName);
			printf("end problem %s\n",problemName);
			printf("\n\n\n"); // nice with blank lines if files are concatenated
			end=1;
			//break;
		}
		else if(!strcmp(keyword,"format")) {
			fscanf(infile,"format version %s",formatVersion);
			printf("format version 1.0\n");
		}
		else if(!strcmp(keyword,"type")) {
			fscanf(infile," = %s",type);
			printf("type = %s\n",type);
			pystr = PyString_FromString(type);
			PyDict_SetItemString(parsedict,"type",pystr);
		}
		else if(!strcmp(keyword,"date")) {
			fscanf(infile," = %s",date);
			printf("date = %s\n",date);
			pystr = PyString_FromString(date);
			PyDict_SetItemString(parsedict,"date",pystr);
		}
		else if(!strcmp(keyword,"url")) {
			fscanf(infile," = %s",url);
			printf("url = %s\n",url);
			pystr = PyString_FromString(url);
			PyDict_SetItemString(parsedict,"url",pystr);
		}
		else if(!strcmp(keyword,"date")) {
			fscanf(infile," = %s",date);
			printf("date = %s\n",date);
		}
		else if(!strcmp(keyword,"reaction")) {   // borde detta vara med alls?
			fscanf(infile," has name = %s",rName);
			printf("reaction_%d\n",keyIndex);
			printf("has name = %s\n",rName);
			if(!strcmp(rName,"uniMolecularMassAction")) {
				fscanf(infile," has localVariableName_%d = %s",&lvIndex,lvName);
				fscanf(infile," has localParameterName_%d = %s",&lpIndex,lpName);
				fscanf(infile," has equation = %s",equation);
				printf("has localVariableName_%d = %s\n",lvIndex,lvName);
				printf("has localParameterName_%d = %s\n",lpIndex,lpName);
				printf("has equation = %s\n",equation);
			}
			else if(!strcmp(rName,"biMolecularMassAction")) {
				fscanf(infile," has localVariableName_%d = %s",&lvIndex,lvName);
				fscanf(infile," has localVariableName_%d = %s",&lvIndex2,lvName2);
				fscanf(infile," has localParameterName_%d = %s",&lpIndex,lpName);
				fscanf(infile," has equation = %s",equation);
				printf("has localVariableName_%d = %s\n",lvIndex,lvName);
				printf("has localVariableName_%d = %s\n",lvIndex2,lvName2);
				printf("has localParameterName_%d = %s\n",lpIndex,lpName);
				printf("has equation = %s\n",equation);
			}
		}
		else if(!strcmp(keyword,"variable")) {
			fscanf(infile," has name = %s is %s",name,property);
			printf("variable_%d has name = %s is %s\n",keyIndex,name,property);
			pydict = PyDict_New();
			pystr = PyString_FromString(name);
			PyDict_SetItemString(pydict,"name",pystr);
			pystr = PyString_FromString(property);
			PyDict_SetItemString(pydict,"property",pystr);
			sprintf(buf,"%d",keyIndex);
			pystr = PyString_FromString(buf);
			pyint = PyInt_FromString(buf,NULL,10);
			PyDict_SetItemString(pydict,"id",pyint);
			PyList_Append(pyvarlist,pydict);
		}
		else if(!strcmp(keyword,"range")) {
			fscanf(infile," has lowerBound = %lf has upperBound = %lf",&lb,&ub);
			printf("range_%d has lowerBound = %lf has upperBound = %lf\n",keyIndex,lb,ub);
		}
		else if(!strcmp(keyword,"memberOfSet")) {
			printf("memberOfSet_%d",keyIndex);
			while(1) {
				fscanf(infile," variable_%d",&varIndex);
				printf(" variable_%d",varIndex);
				if(endOfLine(infile)) break;
			}
			printf("\n");
		}
		else if(!strcmp(keyword,"possibleReaction")) {
			fscanf(infile," of variable_%d",&varIndex);
			fscanf(infile," has type = %s",reactionType);
			printf("possibleReaction_%d of variable_%d\n",keyIndex,varIndex);
			printf("has type = %s\n",reactionType);
			if(!strcmp(reactionType,"uniMolecularMassAction")) {
				fscanf(infile," has spaceOfVariable X1 = %s",varSetX1);
				fscanf(infile," has rangeOfParameter k1 = %s",rangeK1);
				
				printf("has spaceOfVariable X1 = %s\n",varSetX1);
				printf("has rangeOfParameter k1 = %s\n",rangeK1);
				//printf("\n");
			}
			else if(!strcmp(reactionType,"biMolecularMassAction")) {
				fscanf(infile," has spaceOfVariable X1 = %s",varSetX1);
				fscanf(infile," has spaceOfVariable X2 = %s",varSetX2);
				fscanf(infile," has rangeOfParameter k1 = %s",rangeK1);
				
				printf("has spaceOfVariable X1 = %s\n",varSetX1);
				printf("has spaceOfVariable X2 = %s\n",varSetX2);
				printf("has rangeOfParameter k1 = %s\n",rangeK1);
				//printf("\n");
			}
			//else discard rest of row...
		}
		else if(!strcmp(keyword,"errorFunction")) {
			fscanf(infile," has type = %s",errorType);
			if(!strcmp(errorType,"minusLogLikelihoodPlusLambdaK")) {
				fscanf(infile," has equation = %s",equation);
				fscanf(infile," has lambda = %lf",&lambda);
				printf("errorFunction\n");
				printf("has type = %s\n",errorType);
				printf("has equation = %s\n",equation);
				printf("has lambda = %lf\n",lambda);
				//printf("\n");
				pydict = PyDict_New();
				pystr = PyString_FromString(errorType);
				PyDict_SetItemString(pydict,"type",pystr);
				pystr = PyString_FromString(equation);
				PyDict_SetItemString(pydict,"equation",pystr);
				pyfloat = PyFloat_FromDouble(lambda);
				PyDict_SetItemString(pydict,"lambda",pyfloat);
				PyDict_SetItemString(parsedict,"errorfunc",pydict);
			}
			//else discard rest of row...
		}
		else if(!strcmp(keyword,"experiment")) {
			fscanf(infile," has name = %s has %s",name,object);
			printf("experiment_%d has name = %s has %s\n",keyIndex,name,object);
		}
		else if(!strcmp(keyword,"sample")) {
			fscanf(infile," of experiment_%d",&index);
			fscanf(infile," has time = %lf",&time);
			printf("sample_%d of experiment_%d\n",keyIndex,index);
			printf("has time = %lf\n",time);
			fscanf(infile," has variable_ =");
			printf("has variable_ =");
			while(1) {
				fscanf(infile," %lf",&value);
				printf(" %lf",value);
				if(endOfLine(infile)) break;
			}
			printf("\n");
			fscanf(infile," has sdev of variable_ =");
			printf("has sdev of variable_ =");
			while(1) {
				fscanf(infile," %lf",&value);
				printf(" %lf",value);
				if(endOfLine(infile)) break;
			}
			printf("\n");
		}
		
		if(!strcmp(keyword,"//")) {
			fgets(rowbuffer,999,infile);  // read rest of line
			printf("//%s",rowbuffer);
		}
		else while(getc(infile)!='\n' && !feof(infile)); // always read to eol including the eol
		
		if(endOfLine(infile) && !end) printf("\n");	// pass on a blank line after a statement
		
		/*
		if(!strcmp(keyword,"length")){
			fscanf(infile," of person = %d \n",&l);
			printf("length of person = %d\n",l);
		}
		else if(!strcmp(keyword,"weight")) {
			fscanf(infile," of person = %d \n",&w);
			printf("weight of person = %d\n",w);
		}
		else if(!strcmp(keyword,"person")) {
			fscanf(infile," is %s \n",property);
			printf("person_%d is %s\n",keyIndex,property);
		}
		 */


	}
	fclose(infile);

	//Adding remaining entries 
	PyDict_SetItemString(parsedict,"variables",pyvarlist);
	return parsedict;
}








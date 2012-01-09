/*
* This is a test to make c compatible with python
* Author : Subhodeep Moitra(smoitra@cs.cmu.edu)
* BSD license
*/

#include<Python.h>
#include<stdio.h>

static PyObject* 
spam_system(PyObject *self,PyObject *args) {
  const char *command;
  int sts;

  if(!PyArg_ParseTuple(args,"s",&command))
    return NULL;
  sts  = system(command);
  return Py_BuildValue("i",sts);

}

static PyMethodDef SpamMethods[] = {
    {"system", spam_system, METH_VARARGS,
	"Execute a Shell command"},
    {NULL,NULL,0,NULL}  
};

PyMODINIT_FUNC
initspam(void)
{
  (void) Py_InitModule("spam",SpamMethods);
}


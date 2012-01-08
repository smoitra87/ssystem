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


int main() {

  return 0;
}


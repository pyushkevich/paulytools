/******************************************************************
 * REGISTRY Library                                               *
 ******************************************************************
 * Author:						Paul Yushkevich
 *
 * Date:							1998
 *
 * Description					Text I/O with a key/value/folder organization
 *
 *
 * Dependencies:				None
 ******************************************************************
 * registry.cpp
 *	-----------
 *
 ******************************************************************/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "registry.h"

// Begin namespace
NAMESPACE_PAULY_START

unsigned int Hashtable::hash(const char *s) const{
	// Some good hashing function
	int hash=0;
	while(*s) {
		hash ^= (*s);
		s++;
	}
	
	return hash;
}

int Hashtable::findSpot(const char *key) const {
	int position = hash(key) % capacity;
	int i = 1;
	
	while(data[position] && strcmp(key,data[position]->key)!=0) {
		position = (position+i*i) % capacity;
		i++;
	}
	
	return position;
}

void Hashtable::rehash() {
	Hashentry **oldData = data;
	int oldCap = capacity;
	
	// Reallocate the data
	capacity+=increment;
	data = new Hashentry *[capacity];
	memset(data,0,capacity*sizeof(Hashentry*));
	
	// Put each entry in its place
	for(int i=0;i<oldCap;i++) {
		if(oldData[i]) {
			data[findSpot(oldData[i]->key)] = oldData[i];
		}
	}
	
	delete oldData;
}

Hashtable::Hashtable(int initialCapacity,int capacityIncrement) {
	capacity = initialCapacity;
	increment = capacityIncrement;
	
	data = new Hashentry *[capacity];
	memset(data,0,capacity*sizeof(Hashentry*));
	size = 0;
}

Hashtable::~Hashtable() {
	removeAll();
	delete data;
}

void Hashtable::put(const char *key,const void *value) {
	int index = findSpot(key);
	Hashentry *entry = data[index];
	
	if(entry==NULL) {
		entry = new Hashentry;
		data[index] = entry;
		entry->key = new char[strlen(key)+1];
		strcpy(entry->key,key);
		
		size++;
		if(size > capacity/2)
			rehash();
	}
	entry->value = (void *) value;
}

void *Hashtable::get(const char *key) const {
	int index = findSpot(key);
	Hashentry *entry = data[index];
	return (!entry) ? entry : entry->value;
}

void Hashtable::remove(const char *key) {
	int i = findSpot(key);
	if(data[i]) {
		delete data[i]->key;
		delete data[i];
		data[i] =NULL;
		size--;
	}
}

void Hashtable::destroy(const char *key) {
	int i = findSpot(key);
	if(data[i]) {
		//if(data[i]->value)
		//	delete data[i]->value;
		delete data[i]->key;
		delete data[i];
		data[i] =NULL;
		size--;
	}
}

void Hashtable::destroyAll() {
	for(int i=0;i<capacity;i++) {
		if(data[i]) {
			// if(data[i]->value)
			// 	delete data[i]->value;
			delete data[i]->key;
			delete data[i];
			data[i] =NULL;
		}
	}
	size = 0;
}

void Hashtable::removeAll() {
	for(int i=0;i<capacity;i++) {
		if(data[i]) {
			delete data[i]->key;
			delete data[i];
			data[i] =NULL;
		}
	}
	size = 0;
}

void Hashtable::getKeys(char **returnArray) {
	int idx=0;
	for(int i=0;i<capacity;i++) {
		if(data[i]) {
			returnArray[idx] = data[i]->key;
			idx++;
		}
	}
}

int Hashtable::getSize() {
	return size;
}


int keySortFunc(const void *p1,const void *p2) {
	char *s1 = *((char **)p1);
	char *s2 = *((char **)p2);

	// Compare the two strings up to the brackets - first folders and keys are compared
	// though
	while(1) {
		if(*s1 < *s2)
			return -1;
		if(*s2 < *s1)
			return 1;
		if(*s1==0)
			return 0;

		// Special case - compare past brackets - this routine will even sort things like
		// array[12][14][1]
		if(*s1=='[') {
			char *b1 = strchr(s1,']');
			char *b2 = strchr(s2,']');
			if(b1==NULL || b2==NULL)
				return strcmp(s1,s2);
			*b1 = *b2 = 0;
			int a1 = atoi(s1+1);
			int a2 = atoi(s2+1);
			*b1=*b2=']';

			if(a1 < a2)
				return -1;
			if(a1 > a2)
				return 1;


			s1 = b1;
			s2 = b2;
		}

		s1++;s2++;
	}
}

void Registry::setValue(char *key,const char *value) {
	// See if there is a dot
	char *dot = strchr(key,'.');
	if(dot) {
		// See if there is a folder for pre-dot stuff.
		*dot = 0;
		Registry *f = getFolder(key);

		// Create a folder if one is not present
		if(!f) {
			f = new Registry();
			f->addIfNotFound = addIfNotFound;
			char *fkey = new char[strlen(key)+3];
			strcpy(fkey,"{}");
			strcat(fkey,key);
			table.put(fkey,f);
		}

		*dot = '.';
		f->setStringValue(dot+1,value);
	}
	else {
		char *newKey = new char[strlen(key)+2];
		strcpy(newKey,"=");
		strcat(newKey,key);
		table.put(newKey,value);
	}
}

char *Registry::getValue(char *key) {
	// If folder has dot, delegate to next folder.
	char *dot = strchr(key,'.');
	if(dot) {
		*dot = 0;
		Registry *f = getFolder(key);
		*dot = '.';
		return f ? f->getValue(dot+1) : NULL;
	}
	
	else {
		char *trueKey = new char[strlen(key)+2];
		strcpy(trueKey,"=");
		strcat(trueKey,key);
		char *s = (char *)table.get(trueKey);
		delete trueKey;
		return s;
	}
}

int Registry::getValueKeys(char *keys[]) {
   int keyArrayIdx = 0;

   table.getKeys(keys);
   qsort(keys,table.getSize(),sizeof(char *),keySortFunc);

	for(int i=0;i<table.getSize();i++) {
		if(keys[i][0] == '=') {
			keys[keyArrayIdx] = keys[i]+1;
         keyArrayIdx++;
      }
   }

   keys[keyArrayIdx] = NULL;
   return keyArrayIdx;
}

int Registry::getFolderKeys(char *keys[]) {
   int keyArrayIdx = 0;

   table.getKeys(keys);
   qsort(keys,table.getSize(),sizeof(char *),keySortFunc);

	for(int i=0;i<table.getSize();i++) {
		if(keys[i][0] == '{') {
			keys[keyArrayIdx] = keys[i]+2;
         keyArrayIdx++;
      }
   }

   keys[keyArrayIdx] = NULL;

   return keyArrayIdx;
}

int Registry::getKeyArraySize() {
   int nKeys = table.getSize();
   return nKeys+1;
}


/**
* Write this folder to output stream, at given indent, sorted by key
*/
void Registry::write(FILE *f,int indent) {
	char *ind = new char[indent+1];
	memset(ind,' ',indent);
	ind[indent] = 0;

	char **keys = new char*[table.getSize()];
	table.getKeys(keys);

	// Sort keys using a simple sorting algorithm
	qsort(keys,table.getSize(),sizeof(char *),keySortFunc);

	// Write each key out
	for(int i=0;i<table.getSize();i++) {
		char *key = keys[i];
		// Write out the folder
		if(key[0]=='{') {
			//Move up to actual key
			key+=2;

			// Write fodler name, bracket
			fprintf(f,"%s%s {\n",ind,key);

			// Write out folder contents
			Registry *fld = getFolder(key);
			fld->write(f,indent+3);

			// Write the closing bracket
			fprintf(f,"%s}\n",ind);
		}
		else if(key[0]=='=') {
			// Write a key - value pair
			key++;

			// Write encoded value
			fprintf(f,"%s%s = ",ind,key);
			writeEncodedString(f,getValue(key));
			fprintf(f,";\n");
		}
	}

	delete keys;
}

void Registry::setFlagAddIfNotFound(bool yesno) {
	addIfNotFound = yesno;

	char **keys = new char*[table.getSize()];
	table.getKeys(keys);

	for(int i=0;i<table.getSize();i++) {
		char *key = keys[i];
		if(key[0]=='{') {
			// Move up to actual key
			key+=2;

			// Write out folder contents
			Registry *fld = getFolder(key);
			fld->setFlagAddIfNotFound(yesno);
		}
	}

	delete keys;
}

void Registry::collectKeys(list<string> &keyList,string prefix) {
	char **keys = new char*[table.getSize()];
	table.getKeys(keys);

	// Write each key out
	for(int i=0;i<table.getSize();i++) {
		char *key = keys[i];
		// Write out the folder
		if(key[0]=='{') {
			//Move up to actual key
			key+=2;

			// Write out folder contents
			Registry *fld = getFolder(key);
			fld->collectKeys(keyList,prefix + key + ".");
		}
		else if(key[0]=='=') {
			// Write a key - value pair
			key++;

			// Write encoded value
			keyList.push_back(prefix+key);
		}
	}
	
	delete keys;
}

void Registry::collectAllKeys(list<string> &keyList) {
	collectKeys(keyList,"");
}

// String encoder
void Registry::writeEncodedString(FILE *f,char *s) {
	for(int i=0;i<strlen(s);i++) {
		char c = s[i];
		if(c > ' ' && c <= '~' && c != ';' && c != '\\') {
			fputc(c,f);
		}
		else {
			fprintf(f,"\\%02x",(int)c);
		}
	}
}	

int Registry::hexDigit(char c) {
	if(c >= 'a')
		return 10+c-'a';
	if(c >= 'A')
		return 10+c-'A';
	return c-'0';
}

void Registry::decodeString(char *s) {
	char *dst = s;
	
	while(*s) {
		if(*s == '\\') {
			*dst = hexDigit(*(s+1))*16 + hexDigit(*(s+2));
			s+=3;
		}
		else {
			*dst = *s;
			s++;
		}
		dst++;
	}
	
	*dst = 0;
}

// Read this folder in
char *Registry::read(char *file,ostringstream &oss) {
	while(1) {
		// We are expecting a key - find its beginnning and end
		char *keyStart = file + strspn(file," \t\n\r");
		char *keyEnd = keyStart + strcspn(keyStart," \t\n\r={");
		
		// If key starts with closing brace, we are done
		if(*keyStart=='}' || *keyStart==0)
			return keyStart;
		
		// Create a key string
		char *key = new char[keyEnd-keyStart+1];
		strncpy(key,keyStart,keyEnd-keyStart);
		key[keyEnd-keyStart]=0;
		
		// Check the next character after the key
		char *nextChar = keyEnd + strspn(keyEnd," \t\n\r");
		
		// See if this is a name-value assignment
		if(*nextChar == '=') {
			// We got us a key-value pair
			char *valueStart = nextChar+1+strspn(nextChar+1," \t\n\r");
			char *valueEnd = strchr(valueStart,';');
			if(valueEnd==NULL) {
				oss << "Key " << key << " not followed by semicolon terminated value.\n";
				return NULL;
			}
			
			// Decode the value string
			*valueEnd = 0;
			decodeString(valueStart);
			
			// Set value
			char *value = new char[strlen(valueStart)+1];
			strcpy(value,valueStart);
			
			// Set the value
			setValue(key,value);
			
			// Update the string pointer
			file = valueEnd+1;
		}
		
		else if(*nextChar=='{') {
			// Create a folder
			Registry *sub = new Registry();
			sub->addIfNotFound = addIfNotFound;

			char *folderEnd = sub->read(nextChar+1,oss);
			
			// Folder end should be a closing brace
			if(folderEnd==NULL)
				return NULL;
			
			if(*folderEnd!='}') {
				oss << "Sub-folder " << key << " is terminated by '" << *folderEnd;
				oss << "', not a '}'\n";
				return NULL;
			}
			
			// Need to prepend key with "{}"
			char *newKey = new char[strlen(key)+3];
			strcpy(newKey,"{}");
			strcat(newKey,key);
			delete key;
			
			// Add this folder
			table.put(newKey,sub);
			
			file = folderEnd+1;
		}
	}
}

/**
* Get a named registry folder.
*/
Registry *Registry::getFolder(char *key,bool addIfNull) {
	// If folder has dot, delegate to next folder.
	char *dot = strchr(key,'.');
	if(dot) {
		*dot = 0;
		Registry *f = getFolder(key,addIfNull);
		*dot = '.';
		return f ? f->getFolder(dot+1,addIfNull) : NULL;
	}

	else {
		char *trueKey = new char[strlen(key)+3];
		strcpy(trueKey,"{}");
		strcat(trueKey,key);
		Registry *f = (Registry *)table.get(trueKey);
		if(f==NULL) {
			if(addIfNull) {
				f=new Registry();
				f->addIfNotFound = addIfNotFound;
				table.put(trueKey,f);
			}
			else {
				f = NULL;
			}
		}
		else {
			delete trueKey;
		}
		return f;
	}
}

/**
* See if anything with this key exists
*/
bool Registry::hasKeyPrivate(char *key) {
	// If folder has dot, delegate to next folder.
	char *dot = strchr(key,'.');
	if(dot) {
		*dot = 0;
      bool hasKey = hasKeyPrivate(key);
      if(hasKey) {
         hasKey = getFolder(key)->hasKeyPrivate(dot+1);
      }
     *dot = '.';
     return hasKey;
	}

	else {
		char *trueKey = new char[strlen(key)+3];
		strcpy(trueKey,"{}");
		strcat(trueKey,key);
		if(table.get(trueKey)) {
         delete trueKey;
         return true;
      }

      strcpy(trueKey,"=");
		strcat(trueKey,key);
		if(table.get(trueKey)) {
         delete trueKey;
         return true;
      }

		return false;
	}
}

Registry::Registry() : table(10,40) {
	addIfNotFound = false;
}

Registry::Registry(const char *fname) : table(10,40) {
	readFromFile(fname);
}

Registry::~Registry() {
   table.destroyAll();
}

/**
* Return value is NULL if value not found
*/
const char *Registry::getStringValue(const char *key,const char *defaultValue,...) {
	// Initialize arg list
	va_list val;
	va_start(val,defaultValue);
	vsprintf(tmpKey,key,val);
	va_end(val);
	
	char *value = getValue(tmpKey);
	if(value) 
		return value;
	
	if(addIfNotFound) 
		setValue(tmpKey,defaultValue);

	return defaultValue;
}

bool Registry::cmpStringValue(const char *key,const char *cmpValue,...) {
	// Initialize arg list
	va_list val;
	va_start(val,cmpValue);
	vsprintf(tmpKey,key,val);
	va_end(val);
	
	char *value = getValue(tmpKey);

	if(value == NULL) {
		return (cmpValue==NULL);
	}
	if(cmpValue == NULL) {
		return false;
	}
	return (0==strcmp(value,cmpValue));
}

/**
* Return value is -1 if value not found, 0 if found.
*/
int Registry::getIntValue(const char *key,int defaultValue,...) {
	// Initialize arg list
	va_list val;
	va_start(val,defaultValue);
	vsprintf(tmpKey,key,val);
	va_end(val);
	
	char *value = getValue(tmpKey);
	if(value)
		return atoi(value);
	
	if(addIfNotFound) {
		char tmpVal[24];
		sprintf(tmpVal,"%d",defaultValue);
		setValue(tmpKey,tmpVal);
	}

	return defaultValue;
}

/**
* Return value is -1 if value not found, 0 if found.
*/
double Registry::getDoubleValue(const char *key,double defaultValue,...) {
	// Initialize arg list
	va_list val;
	va_start(val,defaultValue);
	vsprintf(tmpKey,key,val);
	va_end(val);
	
	char *value = getValue(tmpKey);
	if(value) 
		return atof(value);
	
	if(addIfNotFound) {
		char tmpVal[60];
		sprintf(tmpVal,"%lg",defaultValue);
		setValue(tmpKey,tmpVal);
	}

	return defaultValue;
}

/**
* Get a bool value
*/
bool Registry::getBooleanValue(const char *key,bool defaultValue,...) {
	// Initialize arg list
	va_list val;
	va_start(val,defaultValue);
	vsprintf(tmpKey,key,val);
	va_end(val);
	
	char *value = getValue(tmpKey);
	if(value) 
		return (atoi(value) > 0);
	
	if(addIfNotFound) {
		setValue(tmpKey,defaultValue ? "1" : "0");
	}

	return defaultValue;
}

/**
* Place a copy of a string into the registry
*/
void Registry::setStringValue(const char *key,const char *value,...) {
	// Initialize arg list
	va_list val;
	va_start(val,value);
	vsprintf(tmpKey,key,val);
	va_end(val);
	
	// Create new value string 
	char *newValue = new char[strlen(value)+1];
	strcpy(newValue,value);
	setValue(tmpKey,newValue);
}

/**
* Set integer value
*/
void Registry::setIntValue(const char *key,int value,...) {
	// Initialize arg list
	va_list val;
	va_start(val,value);
	vsprintf(tmpKey,key,val);
	va_end(val);
	
	// Create new value string 
	char *str = new char[12];
	sprintf(str,"%d",value);
	setValue(tmpKey,str);
}

/**
* Set double value
*/
void Registry::setDoubleValue(const char *key,double value,...) {
	// Initialize arg list
	va_list val;
	va_start(val,value);
	vsprintf(tmpKey,key,val);
	va_end(val);
	
	// Create new value string 
	char *str = new char[24];
	sprintf(str,"%lg",value);
	setValue(tmpKey,str);
}

/**
* Set boolean value
*/
void Registry::setBooleanValue(const char *key,bool value,...) {
	// Initialize arg list
	va_list val;
	va_start(val,value);
	vsprintf(tmpKey,key,val);
	va_end(val);
	
	// Create new value string 
	char *str = new char[24];
	sprintf(str,"%d",(value) ? 1 : 0);
	setValue(tmpKey,str);
}

Registry &Registry::getSubFolder(const char *key,...) {
	// Initialize arg list
	va_list val;
	va_start(val,key);
	vsprintf(tmpKey,key,val);
	va_end(val);

	return *getFolder(tmpKey);
}

bool Registry::isFolder(const char *key,...) {
	// Initialize arg list
	va_list val;
	va_start(val,key);
	vsprintf(tmpKey,key,val);
	va_end(val);

	return (getFolder(tmpKey,false) != NULL);
}

bool Registry::hasKey(const char *key,...) {
   // Initialize arg list
	va_list val;
	va_start(val,key);
	vsprintf(tmpKey,key,val);
	va_end(val);

	return hasKeyPrivate(tmpKey);
}

/**
* Write the folder to a disk
*/
void Registry::writeToFile(const char *pathname) {
	FILE *f = fopen(pathname,"wt");
	if(!f)
		throw RException("Unable to open file for writing.");
	
	write(f,0);
	fclose(f);
}

/**
* Read folder from file
*/
void Registry::readFromFile(const char *pathname) {
	FILE *f = fopen(pathname,"rb");
	if(!f)
		throw RException("Unable to open file for reading.");
	
	fseek(f,0,SEEK_END);
	int size = ftell(f);
	fseek(f,0,SEEK_SET);

	
	char *buffer = new char[size+1];
	fread(buffer,1,size,f);
	buffer[size] = 0;
	
	ostringstream oss;
	char *rc = read(buffer,oss);
	delete buffer;

	if(rc == NULL)
		throw RException(oss.str());
	
	fclose(f);
}

char *Registry::tmpKey = new char[1024];

/*
void main(void) {
	Registry r;
	
	r.setIntValue("model.figureCount",5);
	for(int i=0;i<5;i++) {
		r.setIntValue("model.figure[%d].primitiveCount",20,i);
		for(int j=0;j<20;j++) {
			r.setDoubleValue("model.figure[%d].primitive[%d].x",((double)rand()) / RAND_MAX,i,j);
			r.setDoubleValue("model.figure[%d].primitive[%d].y",((double)rand()) / RAND_MAX,i,j);
			r.setDoubleValue("model.figure[%d].primitive[%d].s",((double)rand()) / RAND_MAX,i,j);
		}
	}

	r.writeToFile("fodler.dat");
}
*/

// Begin namespace
NAMESPACE_PAULY_END

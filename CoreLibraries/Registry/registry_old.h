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
 * registry.h
 *	-----------
 *
 ******************************************************************/
#ifndef _Registry_H_
#define _Registry_H_

#ifdef _MSC_VER
# pragma warning(disable:4786)  // '255' characters in the debug information
# pragma warning(disable:4244)  // conversion from double to float
# pragma warning(disable:4305)  // conversion from double to float (5.0)
#endif //_MSC_VER

#include <stdio.h>
#include <mylibs.h>
#include <string>

#include <sstream>
#include <iostream>
#include <list>
#include <map>

using namespace std;

// Begin namespace
NAMESPACE_PAULY_START

struct Hashentry {
	char *key;
	void *value;
};

class Hashtable {
private:
	int capacity;
	int increment;
	int size;

	Hashentry **data;

	unsigned int hash(const char *s) const;
	int findSpot(const char *key) const;
	void rehash();

public:	
	Hashtable(int initialCapacity,int capacityIncrement);
	~Hashtable();

	void put(const char *key,const void *value);
	void *get(const char *key) const;
	void remove(const char *key);
	void destroy(const char *key);
	void destroyAll();
	void removeAll();
	void getKeys(char **returnArray);
	int getSize();
};

/******************************************************************
 * Exception associated with registry processing
 ******************************************************************/
class RException {
	string m_error;
public:
	RException(string error) : m_error(error) {}
	RException() {}
	~RException() {}

   void print(ostream &s) const {
      s << m_error;
   }
};

inline ostream& operator << (ostream &s,const RException &exc) {
   exc.print(s);
   return s;
}

/******************************************************************
 * Registry Class
 ******************************************************************/
class Registry {
private:
	Hashtable table;

	static char *tmpKey;

	// A flag as to whether keys and folders that are read and not found should be created and 
	// populated with default values.
	bool addIfNotFound;

	void setValue(char *key,const char *value);
	char *getValue(char *key);
	void write(FILE *f,int indent);
	void writeEncodedString(FILE *f,char *s);
	int hexDigit(char c);
	void decodeString(char *s);
	char *read(char *file,ostringstream &oss);
	Registry *getFolder(char *key,bool addIfNull = true);
    bool hasKeyPrivate(char *key);

	void collectKeys(list <string> &list,string prefix);

public:
	Registry();
	Registry(const char *fname);
   ~Registry();

	const char *getStringValue(const char *key,const char *defaultValue,...);
	int getIntValue(const char *key,int defaultValue,...);
	double getDoubleValue(const char *key,double defaultValue,...);
	bool getBooleanValue(const char *key,bool defaultValue,...);

	bool cmpStringValue(const char *key,const char *cmpValue,...);

	void setFlagAddIfNotFound(bool yesno);

	void setStringValue(const char *key,const char *value,...);
	void setIntValue(const char *key,int value,...);
	void setDoubleValue(const char *key,double value,...);
	void setBooleanValue(const char *key,bool value,...);

	Registry &getSubFolder(const char *key,...);
	bool hasKey(const char *key,...);
	bool isFolder(const char *key,...);

   /**
    * These methods allow us to list the keys in the registry.  First two methods both take a
    * string array that must be allocated to hold getKeyArraySize() strings.  (The caller allocates the array)
    * The first method places all value-bearing keys in the folder into the array and returns the number of keys.
    * The array is NULL-terminated.  
    * The second method returns the array and number of subfolder names in that folder.
    * Keys returned are pointers into the registry data, so they should not be modified/deleted
    */
   int getValueKeys(char *keys[]);
   int getFolderKeys(char *keys[]);
   int getKeyArraySize();

	// Collect keys in a list
	void collectAllKeys(list <string> &list);

   void writeToFile(const char *pathname);
	void readFromFile(const char *pathname);
};

// Begin namespace
NAMESPACE_PAULY_END

#endif

/*
This is the patch-like extension to the original PicoJSON parser, intended to use only by #include at the end of
the picojson/picojson.h 'value' class declaration.

Please do NOT include it separately.
*/

public:

  class Exception : public std::logic_error
  {
    public:
      Exception(const std::string& msg) : std::logic_error(msg) {}
  };

public:
  explicit value(uint64_t i) : type_(int64_type)
  {
        if (i > (uint64_t)std::numeric_limits<int64_t>::max())
          throw Exception("JSON node uint64_t ctor: out of range");

        u_.int64_ = i;
  }

  const value& operator[](const std::string& keyname) const
  {
        if (type_ != object_type)  // is<object>()
          throw Exception(std::string("JSON node operator[]: is not an object. keyname=") + keyname);

        const object& obj = *u_.object_; // get<object>();
        const object::const_iterator it = obj.find(keyname);
        if (it == obj.end())
          throw Exception(std::string("JSON node operator[]: key not found. keyname=") + keyname);

        return it->second;
  }

  value& operator[](const std::string& keyname)
  {
        switch (type_) {
          case null_type:
            *this = value(object());
            //assert(type_ == object_type);
          case object_type:
            break;
          default:
            throw Exception(std::string("JSON node operator[]: is not an object. keyname=") + keyname);
        }

        object& obj = *u_.object_; // get<object>();
        object::iterator it = obj.find(keyname);
        if (it == obj.end())
          return obj[keyname] = picojson::value();  // creates a new empty key

        return it->second;
  }

  const value& operator[](unsigned int index) const
  {
        if (type_ != array_type)  // is<array>()
          throw Exception("JSON node operator[]: is not an array");

        const array& ary = *u_.array_; // get<array>();
        if (index >= ary.size())
          throw Exception("JSON node operator[]: index out of range");

        return ary[index];
  }

  value& operator[](unsigned int index)
  {
        if (type_ != array_type)  // is<array>()
          throw Exception("JSON node operator[]: is not an array");

        array& ary = *u_.array_; // get<array>();
        if (index >= ary.size())
          throw Exception("JSON node operator[]: index out of range");

        return ary[index];
  }

  bool isMember(const std::string& keyname) const
  {
        if (type_ != object_type)  // is<object>()
          throw Exception(std::string("JSON node isMember: is not an object. keyname=") + keyname);

        const object& obj = *u_.object_; // get<object>();
        return obj.find(keyname) != obj.end();
  }

  bool isMember(const char *keyname) const { return isMember(std::string(keyname)); }

  bool isObject() const { return type_ == object_type; }

  bool isArray() const { return type_ == array_type; }

  bool isNumeric() const { return type_ == int64_type || type_ == number_type; }

  bool isIntegral() const { return type_ == int64_type; }

  bool isReal() const { return type_ == number_type; }

  bool isString() const { return type_ == string_type; }

  bool isNull() const { return type_ == null_type; }

  int asInt() const
  {
        if (type_ != int64_type) // is<int64_t>()
          throw Exception("JSON node asInt: has no integral type");

        int64_t value = u_.int64_; // get<int64_t>();
        if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
          throw Exception("JSON node asInt: out of range");

        return static_cast<int>(value);
  }

  unsigned int asUInt() const
  {
        if (type_ != int64_type) // is<int64_t>()
          throw Exception("JSON node asInt: has no integral type");

        int64_t value = u_.int64_; // get<int64_t>();
        if (value < std::numeric_limits<unsigned int>::min() || value > std::numeric_limits<unsigned int>::max())
          throw Exception("JSON node asUInt: out of range");

        return static_cast<unsigned int>(value);
  }

  uint64_t asUInt64() const
  {
        if (type_ != int64_type) // is<int64_t>()
          throw Exception("JSON node asInt: has no integral type");

        return static_cast<uint64_t>(u_.int64_); // get<int64_t>();
  }

  double asDouble() const
  {
        switch (type_) {
          case number_type:
            return u_.number_;
          case int64_type:
            return u_.int64_;
          default:
            throw Exception("JSON node asDouble: has no numeric type");
        }
  }

  float asFloat() const
  {
        switch (type_) {
          case number_type:
            return u_.number_;
          case int64_type:
            return u_.int64_;
          default:
            throw Exception("JSON node asFloat: has no numeric type");
        }
  }

  const std::string& asString() const
  {
        if (type_ != string_type) // is<std::string>()
          throw Exception("JSON node asString: has not a string type");

        return *u_.string_; // get<std::string>();
  }

  uint64_t asBigUInt() const
  {
    if (isIntegral())
        return asUInt64();

    if (isReal())
        return std::max(ceil(asDouble()),0.0);

    if (isString()) {
        if (asString().empty() || asString()[0] == '-')
            throw Exception("big number recieved negative or no value");
        return strtoul(asString().c_str(), NULL, 0);
  }
    throw Exception("big number string not interpretable");
  }

  bool asBool() const
  {
        if (type_ != boolean_type) // is<int64_t>()
          throw Exception("JSON node asBool: has no boolean type");

        return u_.boolean_; // get<bool>();
  }

  const array& asArray() const
  {
        if (type_ != array_type)  // is<array>()
          throw Exception("JSON node asArray: is not an array");

        return *u_.array_; // get<array>();
  }

  array& asArray()
  {
        if (type_ != array_type)  // is<array>()
          throw Exception("JSON node asArray: is not an array");

        return *u_.array_; // get<array>();
  }

  const object& asObject() const
  {
        if (type_ != object_type)  // is<object>()
          throw Exception("JSON node asObject: is not an object");

        return *u_.object_; // get<object>();
  }

  object& asObject()
  {
        if (type_ != object_type)  // is<object>()
          throw Exception("JSON node asObject: is not an object");

        return *u_.object_; // get<object>();
  }

  value& operator=(int v) { return *this=value((int64_t)v); }
  value& operator=(unsigned int v) { return *this=value((uint64_t)v); }
  value& operator=(const char *v) { return *this=value(std::string(v)); }

  #define OP_ASSIGN(type) \
  value& operator=(type v) { return *this=value((type)v); }

  OP_ASSIGN(const std::string&);
  OP_ASSIGN(uint64_t);
  OP_ASSIGN(float);
  OP_ASSIGN(double);
  OP_ASSIGN(bool);
  OP_ASSIGN(const array&);
  OP_ASSIGN(const object&);

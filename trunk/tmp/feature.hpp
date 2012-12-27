#ifndef FEATURE_HPP
#define FEATURE_HPP

class Feature {

public:

  Feature();
  ~Feature();

  bool isNumerical() const;
  bool isCategorical() const;
  bool isTextual() const;

protected:

  enum Type { NUM, CAT, TXT, UNKNOWN };

  virtual initialize() = 0;

private:

  Type type_;

};

class NumFeature : public Feature {

public:

  NumFeature();
  ~NumFeature();

protected:

  virtual initialize();

private:



};


#endif

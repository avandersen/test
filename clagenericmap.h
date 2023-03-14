#pragma once
#include "ecmiscbasic.h"

struct DEntry;

class ClAGenericMap {
public:
  ECVARLISTEMPTY(ClAGenericMap)
    ECSUBLIST(ClAGenericMap, (ClADataMap)(ClDepthCorrection))
    ECSERIALIZE(ClAGenericMap)
    /*  virtual size_t remove(const double &key) =0;
    virtual QList<double> keys() const = 0;
    QList<double> uniqueKeys() const = 0;
    QList<boost::any> values() const = 0;
    QList<boost::any> values(const double &key) = 0;
    const boost::any value(const double &key, const boost::any& def = boost::any()) const = 0;
    size_t count() const = 0;
    size_t count(const double& key) const = 0;
    bool contains(const double& key) const = 0;*/
  EcUuid _dbid = EcUuid::sequential();

    static ClAGenericMap_PT createMap(ClDataType dt);

  virtual ~ClAGenericMap() {}
  virtual bool isInterval() const { return false; }
  virtual size_t count() const = 0;
  virtual size_t count(const double& key) const = 0;
  virtual bool contains(const double& key) const = 0;
  virtual bool empty() const = 0;
  virtual void clear() = 0;
  virtual ClDataType dataType() const = 0;
  virtual QList<double> keys() const = 0;
  virtual QList<double> uniqueKeys() const = 0;
  virtual QList<boost::any> values() const = 0;
  virtual QList<boost::any> values(const double& key) const = 0;
  virtual QList<DEntry> entries(const double& key) const = 0;
  virtual const boost::any value(const double& key, const boost::any& def = boost::any()) const = 0;
  virtual ClAGenericMap_PT getCopy() const = 0;
  virtual ClAGenericMap_PT getCopyNewIds() const = 0;
  virtual ClADataMap_PT asDataMap() = 0;

  virtual double getTopDepth() const = 0;
  virtual double getBottomDepth() const = 0;
  virtual double getBottomExtent() const { return getBottomDepth(); }
  virtual double getTopExtent() const { return getTopDepth(); }
};
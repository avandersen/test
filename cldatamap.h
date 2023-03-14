#pragma once

#include "clagenericmap.h"
#include "clsql.h"
#include "ecdebug.h"
#include "ecmultiindexcontainer.h"
#include "IntervalTree.h"

double capPrecision(double v);

// interval no overlap (any)
// interval overlap (composite, must have bottomdepth)
// single depth, single value (any)
// single depth, multiple values (any except empty)
// discontinuous (double)

class ClCValue {
public:

  ECVARLIST(ClCValue,
    ((boost::variant<boost::blank, QMap<QString, boost::any>, int, double, QString>), _vmap, {}))
    //    ((EcVariant), _vmap, {}))
    ECSERIALIZE(ClCValue)


  ClCValue();
  ClCValue(int v);
  ClCValue(double v);
  ClCValue(QString v);

  struct ClCValue_iterator {
    bool _isfirst = true;
    const ClCValue* _me;
    QMap<QString, boost::any>::const_iterator _it;
    bool _isend = false;
    ClCValue_iterator& operator++();
    bool operator==(const ClCValue_iterator &) const;
    bool operator!=(const ClCValue_iterator &) const;
    const QString &key();
    boost::any valueAny();

  };
  void dumpContents();

  template<class T> bool isBaseValue() const;
  bool hasValueAny(const QString& str) const;
  template<class T> bool hasValue(const QString& str) const;
  template<class T> bool hasValueBase() const;
  boost::any valueAny(const QString& str) const;
  template<class T> T value(const QString& str) const;
  template<class T> T valueBase() const;
  void setValueAny(const QString& str, const boost::any& v);
  template<class T> void setValue(const QString& str, const T& v);
  void clearValue(const QString& str);
  QMap<QString, boost::any> asMap() const; // AVA: legacy use builtin iterator
  bool operator ==(const ClCValue& other) const;
  bool operator !=(const ClCValue& other) const;
  ClCValue& operator=(int v);
  ClCValue& operator=(double v);
  ClCValue& operator=(const QString& v);
  int count() const;
  boost::any getIndexedAny(int i);
  bool hasPointSet() const;
  bool hasPosition() const;
  double pointSetTopDepth() const;
  double pointSetBottomDepth() const;
  void shiftDepths(double delta);

  ClCValue_iterator begin() const;
  ClCValue_iterator end() const;
  
  template<class T> void dbgCheckType(const QString&) const;

/*  bool hasValueAny(const QString& str) const { return _vmap.contains(str); }
  template<class T> bool hasValue(const QString &str) const { return _vmap.value(str).type() == typeid(T); }
  boost::any valueAny(const QString &str) const { return _vmap.value(str); }
  template<class T> T value(const QString &str) const { return anyc<T>(_vmap.value(str)); }
  void setValueAny(const QString &str, const boost::any &v) { _vmap[str]=v; }
  template<class T> void setValue(const QString &str, const T &v) { _vmap[str] = cr_any<T>(v); }
  void clearValue(const QString& str) { _vmap.remove(str); }
  QMap<QString, boost::any> asMap() const;
  bool operator ==(const ClCValue& other) const;
  bool operator !=(const ClCValue& other) const { return !(*this == other); }*/
};


struct DEntry {
  DEntry() {}
//  DEntry(const std::pair<double, boost::any>& v);
  DEntry(double top, double bottom, const boost::any&, const EcUuid&);
  DEntry(double top, const boost::any&);
  ECVARLIST(DEntry,
    (double, first, 0)
    (double, bottom, 0)
    (EcUuid, dmid, {})
    (EcUuid, uuid, {})
    (boost::any, second, {}))
  ECSERIALIZE(DEntry)
  bool operator<(const DEntry& other) const { return uuid < other.uuid; }
  int interval_group() const {
    return (first - bottom == 0) ? 0 : MAX(1, MIN(7, ceil(log10(fabs(first - bottom)))));
  }

  template<bool isread, class V> struct serialize;
  void serializeEnabled() {}
};

#if 1

#define USEMULTIINDEXDM
ECMIC(ClDataMapDataContainer, DEntry,
(boost::multi_index_container<DEntry, 
  bmix::indexed_by<
    bmix::ordered_non_unique<bmix::tag<tag_depth>, bmix::member<DEntry, double, & DEntry::first> >,
    bmix::ordered_non_unique<bmix::tag<tag_ivgroup>,
      bmix::composite_key<DEntry,
        bmix::const_mem_fun<DEntry, int, &DEntry::interval_group>,
        bmix::member<DEntry, double, & DEntry::first>
      >
    >,
    bmix::ordered_unique<bmix::tag<tag_id>, bmix::member<DEntry, EcUuid, &DEntry::uuid> >,
    bmix::ordered_non_unique<bmix::tag<tag_bottom>, bmix::member<DEntry, double, &DEntry::bottom> >
  >
 >), (2)(2)(2)(2))

struct ClDataMapRange;

struct ClDataMapRangeIterator {
  ClDataMapRange* dr = NULL;
  int curix = 0;
  decltype(ClDataMapDataContainer().get<tag_ivgroup>().begin()) cur;
  ClDataMapRangeIterator& operator++();
  bool operator==(const ClDataMapRangeIterator& other);
  bool operator!=(const ClDataMapRangeIterator& other);
  const DEntry& operator*() const;
  const DEntry* operator->() const;
  ClDataMapRangeIterator() {};
};

struct ClDataMapRange {
  const ClADataMap* dm;
  double top, bot;
  std::vector<std::pair<decltype(ClDataMapDataContainer().get<tag_ivgroup>().begin()), decltype(ClDataMapDataContainer().get<tag_ivgroup>().begin())> > range;
  ClDataMapRangeIterator begin() {
    ClDataMapRangeIterator res;
    if (range.size() == 0) {
      res.curix = -1;
    }
    else {
      res.curix = 0;
      res.cur = range[0].first;
    }
    res.dr = this;
    return res;
  }
  ClDataMapRangeIterator end() {
    ClDataMapRangeIterator res;
    res.curix = -1;
    res.dr = this;
    return res;
  }
};


//typedef ClDBDataMap ClDataMapDataContainer;
#else
typedef std::multimap<double, boost::any> ClDataMapDataContainer;
#endif
class ClADataMap : public ClAGenericMap {
  //  boost::any &	operator[] ( const double &key ) { return QMap<double, boost::any>::operator [](key); }
  //  const boost::any	operator[] ( const double &key ) const  { return QMap<double, boost::any>::operator [](key); }
public:

  ECVARLIST(ClADataMap,
    (ClDataMapDataContainer, _data, ClDataMapDataContainer())
    (QString, _unitName, {})
    (ClDataType, _dataType, CDT_None))

    ECSUPERLIST(ClADataMap, (ClAGenericMap))
    ECSERIALIZE(ClADataMap)
    //  double _bulkShift;
    ClADataMap() { _dataType = CDT_None; } // AVA: TBD Remove

//  EcIntervalTree<double, boost::any> _intervalTree = { 16, -100000000, 1000000000 };


  typedef decltype(_data.get<tag_depth>().begin()) iterator;
  typedef decltype(_data.get<tag_depth>().cbegin()) const_iterator;
  iterator getFirstIteratorAbove(double);
  const_iterator getFirstIteratorAbove(double) const;
  iterator begin() { return _data.get<tag_depth>().begin(); }
  iterator end() { return _data.get<tag_depth>().end(); }
  const_iterator begin() const { return _data.get<tag_depth>().cbegin(); }
  const_iterator end() const { return _data.get<tag_depth>().cend(); }

  iterator erase(iterator it);
  virtual size_t remove(const double &key);
#ifdef WIN32
  // only here to prevent implicit casts
  template<class T> iterator insert(const double &key, const T& data) { static_assert(!std::is_same<T,T>::value); }
#endif
  void insertAtEnd(const double& key, boost::any&& data);
  iterator insert(const double& key, boost::any&& data);
  iterator insert(const double &key, const boost::any& data);
  iterator insert(const DEntry &entry);
#ifdef WIN32
  // only here to prevent implicit casts
  template<class T> iterator insertMulti(const double &key, const T& data) { static_assert(!std::is_same<T, T>::value); }
#endif
  ClDataMapRange getRange(double top, double bottom) const;
  std::vector<std::pair<double, double> > getDataGapsInInterval(double top, double bot) const;
  void insertDEntryList(const QList<DEntry>& data);
  size_t insertMulti(const double &top, const double &bot, const boost::any& data, int ix);
  size_t insertMulti(const double &top, const double &bot, const boost::any& data);
  void insertMultiList(const double &top, const double &bot, const QList<boost::any>& data);
  iterator replace(const double& top, const double & bottom, const boost::any& value);
  iterator replace(const double& top, const boost::any& value); // AVATBD REMOVE
  iterator replace(iterator& it, boost::any val) {
#ifdef USEMULTIINDEXDM
    auto tmp = *it;
    tmp.second = val;
    _data.get<tag_depth>().replace(it, tmp);
#else
    it->second = val;
#endif
    return it;
  };
#ifdef WIN32
  // only here to prevent implicit casts to any... use cr_any<tp>(val)
  template<class T>
  void setValueUpwards(const double &key, const T& value) { static_assert(!std::is_same<T, T>::value); }
  template<class T>
  void setValueDownwards(const double &key, const T& value) { static_assert(!std::is_same<T, T>::value); }
  template<class T>
  void setSingleValue(const double &key, const boost::any & value) { static_assert(!std::is_same<T, T>::value); }
#endif
  std::vector<boost::any> setSingleValue(const double &key, const boost::any & value);
  std::vector<boost::any> setValueUpwards(const double &key, const boost::any& value);
//  ClDataModifications ClADataMap::setValueUpwardsDM(const double &key, const boost::any& value) {
  std::vector<boost::any> setValueDownwards(const double &key, const boost::any& value);
  std::vector<boost::any> setValueDir(const double &key, const boost::any& value, bool upwards);
  std::vector<boost::any> replaceValues(const double &key, const std::vector<boost::any>&);

  std::set<EcUuid> uuidsAtDepth(double depth) const;
  boost::any getValueUpwards(const double &key) const;
  boost::any getValueDownwards(const double &key)const;
  boost::any getValueDir(const double &key, bool upwards)const;
  double getContentsWidth(ClDocument_CPT cldoc, ClIImageColumn_PT icol);

  iterator find(const double& key);
  iterator findMulti(const double& key, int ix);
  const_iterator find(const double& key) const;
  //  const_iterator constFind(const double &key) const { return _data.constFind(key); }
  iterator lowerBound(const double &key);
  const_iterator lowerBound(const double &key) const;
  iterator upperBound(const double &key);
  const_iterator upperBound(const double &key) const;
  int getIteratorIndex(const_iterator& it);

  const boost::any valueApproximate(const double &key, const boost::any& def = boost::any()) const;
  QList<boost::any> valuesApproximate(const double &key) const;
  virtual size_t countApproximate(const double& key) const;
  virtual ClADataMap_PT extractInterval(double top, double bottom) const; // works with bidir double point data
  bool isSuperSetOf(ClADataMap_PT);

  ClADataMap(const ClDepthCorrection &);
  ClADataMap(ClDataType, const EcUuid& id);
  ClADataMap(ClDataType, std::set<std::pair<double, int> >, boost::any, const EcUuid &id);
  virtual ~ClADataMap() {}

  bool isDefinedAt(double depth) const;
  virtual ClADataMap_PT asDataMap() { return dPTc<ClADataMap>(getCopy()); }

  // AVA: only for CDT_Depth_Map
  //  double bulkShift() const { return _bulkShift; } 
  //  void setBulkShift(double bs) { _bulkShift=bs; }

  double getFirstEmptyDepthBelow(double dpth) const; // interval data
  double getFirstEmptyDepthAbove(double dpth) const;

  // AVA: assumes existence
  double getMaxValue() const;
  double getMinValue() const;
  void getValueRange(double& absmin, double& absmax, double& rngmin, double& rngmax) const;
  double getDouble(double depth) const;
  int getInt(double depth) const;
  QColor getColor(double depth) const;
  QString unitName() const { return _unitName; }
  void setUnitName(const QString &str) { _unitName = str; }
  bool isIntervalOverlap(double top, double bottom) const;

  bool isIntervalEmpty(double top, double bottom) const;

  // AVA: always includes top and bottom
  std::set<double> getKeysInInterval(double top, double bottom, bool includefirstoutside) const;
  ClADataMap_PT getBidirValuesInInterval(double top, double bottom);
  double getDualDataPosAbove(double depth, bool *valid, bool exclusive = false) const;
  double getDualDataPosBelow(double depth, bool *valid, bool exclusive = false) const;
  double getDataPosAboveBottom(double vp, bool* valid, bool exclusive= false) const;
  double getDataPosBottom(double vp, bool* valid) const;
  double getDataPosAbove(double depth, bool *valid, bool exclusive = false) const;
  double getDataPosBelow(double depth, bool *valid, bool exclusive = false) const;
  double getClosestIntervalBorder(double vp, bool* ok=NULL) const;
  double getClosestDataPos(double depth, bool *ok=NULL) const;
  double getDoubleInterpolated(double depth, bool *ok = NULL, bool preferdownward=true) const;
//  void getTightDataRange(double top, double bottom, double* maxtop, double* maxbottom);

  void rwXml(const QString &key, QDomElement &, bool isRead, bool isArchive, bool ignoreimages, const QMap<QString, QString> &imagemapnames);
  void purgeImageData();
  void reloadImageData();

  // ClAGenericMap methods below
  virtual bool isInterval() const { return (_dataType == CDT_Age_Interval || _dataType == CDT_Int_Interval || _dataType == CDT_Color_Interval || _dataType == CDT_String_Interval || _dataType == CDT_Lithology_Interval); }
  virtual size_t count() const { return _data.size(); }
  virtual size_t count(const double& key) const { return _data.get<tag_depth>().count(capPrecision(key)); }
  virtual bool contains(const double& key) const { return _data.get<tag_depth>().find(capPrecision(key))!=_data.get<tag_depth>().end(); }
  virtual bool containsIntervalBorder(const double& key) const { return _data.get<tag_depth>().find(capPrecision(key)) != _data.get<tag_depth>().end() || _data.get<tag_bottom>().find(capPrecision(key)) != _data.get<tag_bottom>().end(); }
  virtual bool empty() const { return _data.empty(); }
  virtual void clear();
  virtual ClDataType dataType() const { return _dataType; }
  virtual QList<double> keys() const;
  QList<double> uniqueKeys() const;
  virtual QList<boost::any> values() const;
  virtual QList<boost::any> values(const double &key) const;
  virtual QList<DEntry> entries(const double& key) const;
  virtual std::set<DEntry> entriesSet(const double& key) const;
  const boost::any value(const double &key, const boost::any& def = boost::any()) const;
  virtual ClAGenericMap_PT getCopy() const { return ClAGenericMap_PT(new ClADataMap(*this)); }
  virtual ClAGenericMap_PT getCopyNewIds() const;
  virtual double getTopDepth() const;
  virtual double getBottomDepth() const;
  virtual double getBottomExtent() const ;
  virtual double getTopExtent() const;

  template<bool isread, class V> struct serialize;
  void serializeEnabled() {}

  friend class ClDataSet;
};

class ClIntervalMap;

class ClIntervalMapIterator {
  const ClIntervalMap *_ivmap;
  ClADataMap::iterator dmit;
  ClIntervalMapIterator(const ClIntervalMap *ivmap, ClADataMap::iterator it) : _ivmap(ivmap), dmit(it) {};
public:
  double topDepth() const;
  double bottomDepth() const;
  boost::any value() const;
  boost::any namedProperty(const QString &) const;

  ClIntervalMapIterator& operator++();
  ClIntervalMapIterator& operator--();
  bool operator==(const ClIntervalMapIterator& other) const;
  bool operator!=(const ClIntervalMapIterator& other) const;

  friend class ClIntervalMap;
};

class ClIntervalMap {
  ClIntervalMapIterator iterator;

  ClADataMap *_dm;
public:
  ClIntervalMapIterator begin() const;
  ClIntervalMapIterator end() const;

  friend class ClIntervalMapIterator;
};

class ClPointsMap : public ClADataMap {
public:
};


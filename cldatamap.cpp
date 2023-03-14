#include "clpch.h"

#include <algorithm>

#include "cldata.h"
#include "cldatamap.h"
#include "cldepthcorrection.h"
#include "cliimagecolumn.h"
#include "climage.h"
#include "clinotescolumn.h"
#include "cldoccmd.h"

#include "eciimage.h"

//DEntry::DEntry(const std::pair<double, boost::any>& v) : first(v.first), bottom(v.first), second(v.second), uuid(EcUuid::sequential()) {
//}
DEntry::DEntry(double t, double b, const boost::any& v, const EcUuid& id) : first(t), bottom(b),  uuid(id), second(v) {
}
DEntry::DEntry(double t, const boost::any& v) : first(t), bottom(t), uuid(EcUuid::sequential()), second(v) {
}

ClDataMapRangeIterator& ClDataMapRangeIterator::operator++() {
  do {
    cur++;
    if (cur == dr->range[curix].second) {
      ++curix;
      if (curix >= (int) dr->range.size()) {
        curix = -1;
      }
      else {
        cur = dr->range[curix].first;
      }
    }
  } while (curix != -1 && (cur->first > dr->bot || cur->bottom < dr->top));
  return *this;
}
bool ClDataMapRangeIterator::operator==(const ClDataMapRangeIterator& other) {
  if (curix != other.curix) return false;
  if (curix == -1) return true;
  return cur == other.cur;
}
bool ClDataMapRangeIterator::operator!=(const ClDataMapRangeIterator& other) { return !(*this == other); }
const DEntry& ClDataMapRangeIterator::operator*() const { return *cur; }
const DEntry* ClDataMapRangeIterator::operator->() const { return &*cur; }

ClDataMapRange ClADataMap::getRange(double top, double bot) const {
  ClDataMapRange res;
  res.dm = this;
  res.top = top;
  res.bot = bot;
  for (int i = 0; i <=7; ++i) {
    double diff = pow(10, i);
    auto rng = _data.get<tag_ivgroup>().range(boost::make_tuple(i, top - diff) <= boost::lambda::_1, boost::make_tuple(i, bot) >= boost::lambda::_1);
    while (rng.first != rng.second && (rng.first->bottom< top || rng.first->first >bot)) rng.first++;
    if (rng.first != rng.second) res.range.push_back(rng);
  }
  return res;
}

ClADataMap::ClADataMap(const ClDepthCorrection& dc) {
  _dbid = dc._dbid;
  _dataType = CDT_Depth_Map;
  //  _bulkShift=dc._bulkshift;
  for (auto it = dc._correctionMap.get<tag_left>().begin(); it != dc._correctionMap.get<tag_left>().end(); ++it) {
    insertMulti(it->left, it->left,boost::any(it->right));
  }
}

std::set<EcUuid> ClADataMap::uuidsAtDepth(double depth) const {
  std::set<EcUuid> res;
  for (auto it = find(depth); it != end() && it->first==depth; ++it) {
    res.insert(it->uuid);
  }
  return res;
}

ClADataMap::ClADataMap(ClDataType cdt, const EcUuid& id) {
  _dbid = id;
  _dataType = cdt;
  //  _bulkShift=0;
}

ClADataMap::ClADataMap(ClDataType cdt, std::set<std::pair<double,int > > sel, boost::any v, const EcUuid& id) {
  _dbid = id;
  for (std::set<std::pair<double, int> >::iterator it = sel.begin(); it != sel.end(); ++it) {
    replace(it->first, v);
  }
  _dataType = cdt;
  //  _bulkShift=0;
}

bool ClADataMap::isDefinedAt(double depthp) const {
  auto rng=getRange(depthp, depthp);

  return rng.begin() != rng.end();
}

int ClADataMap::getIteratorIndex(const_iterator& it) {
  if (it == end()) return 0;
  auto ti = find(it->first);
  int ix = 0;
  while (true) {
    if (ti == end() || ti->first != it->first) return 0;
    if (ti == it) return ix;
    ++ix;
    ++ti;
  }
  return 0;
}

size_t ClADataMap::countApproximate(const double& key) const {
  size_t ix = 0;
  for (auto it = _data.get<tag_depth>().lower_bound(key-0.0000001); it != end() && it->first < key+0.0000001; ++it) {
    ++ix;
  }
  return ix;
}

const boost::any ClADataMap::valueApproximate(const double &key, const boost::any& def) const {
  size_t ix = 0;
  if (contains(key)) return value(key);
  auto it = _data.get<tag_depth>().lower_bound(key-0.0000001);
  for (; it != end() && it->first < key+0.0000001; ++it) {
    ++ix;
    break;
  }

  if (ix) return it->second;
  return def;
}

QList<boost::any> ClADataMap::valuesApproximate(const double &key) const {
  QList<boost::any> res;
  auto it = _data.get<tag_depth>().lower_bound(key-0.0000001);
  for (; it != end() && it->first < key+0.0000001; ++it) {
    res.push_back(it->second);
  }

  return res;
}
double ClADataMap::getContentsWidth(ClDocument_CPT cldoc, ClIImageColumn_PT icol) {
  double res = 1;
  for (ClADataMap::iterator it = begin(); it != end(); ++it) {
    EcIImage_PT im = icol->getPositionedImage(it, cldoc, it->first, icol->_stretchToPreserveWidthOnIrregularShifts, icol->_fixedHorizontalScale);
    if (im) res = MAX(res, im->visualBoundingRect().width());
    im = icol->getPositionedImage(it, cldoc, it->bottom, icol->_stretchToPreserveWidthOnIrregularShifts, icol->_fixedHorizontalScale);
    if (im) res = MAX(res, im->visualBoundingRect().width());
  }
  return res;
}

bool ClADataMap::isIntervalOverlap(double topp, double bottomp) const {
  auto rng = getRange(topp, bottomp);
  for (auto it = rng.begin(); it != rng.end(); ++it) {
    if (!(topp >= it->bottom || bottomp <= it->first)) return true;
  }
  return false;
}

/*  double top = capPrecision(topp);
  double bottom = capPrecision(bottomp);
  bool atopvalid,btopvalid;
  double abovetop = getDataPosAbove(top, &atopvalid, false);
  getDataPosBelow(top, &btopvalid, true);
  if (atopvalid && !value(abovetop).empty() && btopvalid) return true;
  bool abottomvalid;
  double abovebottom = getDataPosAbove(bottom, &abottomvalid, true);
  if (dataType() == CDT_Image) {
    if (abottomvalid && (abovebottom > top || (!value(abovebottom).empty() && find(abovebottom)->bottom > top))) return true;
  }
  else {
    if (abottomvalid && abovebottom > top) return true;
  }
  return false;
}*/


bool ClADataMap::isIntervalEmpty(double topp, double bottomp) const {
  double top = capPrecision(topp);
  double bottom = capPrecision(bottomp);

  auto rng = getRange(top, bottom);
  for (auto r : rng) {
    if (!(r.first >= bottomp || r.bottom <= top)) return false;
  }
  return true;
}

void ClADataMap::getValueRange(double &absmin, double &absmax, double &rngmin, double &rngmax) const {
  // SJM 2023-01-13: this function replaces getMinValue and getMaxValue, and also finds a good range for 2d gradient
  bool found = false;
  double curmin = 0;
  double curmax = 0;
  for (const_iterator it = begin(); it != end(); ++it) {
    if (it->second.type() == typeid(double)) {
      double v = anyc<double>(it->second);
      if (!found) {
        found = true;
        curmin = curmax = v;
      }
      else if (v < curmin) {
        curmin = v;
      }
      else if (v > curmax) {
        curmax = v;
      }
    }
    if (it->second.type() == typeid(ClCValue) && anyc<ClCValue>(it->second).hasValue<std::vector<double>>(VALUE)) {
      std::vector<double> dv = anyc<ClCValue>(it->second).value<std::vector<double>>(VALUE);
      if (dv.size() == 0) continue;
      // SJM 2023-01-06: exclude undefined values from check
      for (auto v : dv) {
        if (v != -999.25) {
          if (!found) {
            curmin = curmax = v;
            found = true;
          }
          else if (v < curmin) {
            curmin = v;
          }
          else if (v > curmax) {
            curmax = v;
          }
        }
      }
    }
  }
  absmin = curmin;
  absmax = curmax;

  if (dataType() == CDT_Double_Vector_Point) {
    size_t rngcnt[50] = { 0 };
    size_t totcnt = 0;
    double w = absmax - absmin;
    if (w == 0) {
      return;
    }
    for (const_iterator it = begin(); it != end(); ++it) {
      if (it->second.type() == typeid(ClCValue) && anyc<ClCValue>(it->second).hasValue<std::vector<double>>(VALUE)) {
        std::vector<double> dv = anyc<ClCValue>(it->second).value<std::vector<double>>(VALUE);
        if (dv.size() == 0) continue;
        // SJM 2023-01-06: exclude undefined values from check
        for (auto v : dv) {
          if (v != -999.25) {
            int n = (int)(50.0 * (v - absmin) / w);
            if (n > 49) n = 49;
            ++rngcnt[n];
            ++totcnt;
          }
        }
      }
    }
    double ratio = 0.0;
    bool fmin = false;
    bool fmax = false;
    for (size_t i = 0; i < 50; i++) {
      ratio += (double)rngcnt[i] / (double)totcnt;
      if (!fmin && ratio > 0.05) {
        fmin = true;
        if (i > 3 && i < 46) rngmin = i * 0.02;
      }
      if (!fmax && ratio > 0.95) {
        fmax = true;
        if (i > 3 && i < 46) rngmax = i * 0.02;
      }
    }
  }
}

double ClADataMap::getMinValue() const {
  bool found = false;
  double curval = 0;
  for (const_iterator it = begin();it != end(); ++it) {
    if (it->second.type() == typeid(double)) {
      if (!found) {
        found = true;
        curval = anyc<double>(it->second);
      }
      else if (curval > anyc<double>(it->second)) {
        curval = anyc<double>(it->second);
      }
    }
    if (it->second.type() == typeid(ClCValue) && anyc<ClCValue>(it->second).hasValue<std::vector<double>>(VALUE)) {
      std::vector<double> dv = anyc<ClCValue>(it->second).value<std::vector<double>>(VALUE);
      if (dv.size() == 0) continue;
      // SJM 2023-01-06: exclude undefined values from check
      for (auto v : dv) {
        if (v != -999.25) {
          if (!found || v < curval) {
            curval = v;
            found = true;
          }
        }
      }
/*
      if (!found) {
        std::vector<double> dv = anyc<ClCValue>(it->second).value<std::vector<double>>(VALUE);
        if (dv.size() == 0)continue;
        found = true;
        if (dv.size()>0) curval = *std::min_element(dv.begin(), dv.end());
      }
      else {
        std::vector<double> dv = anyc<ClCValue>(it->second).value<std::vector<double>>(VALUE);
        if (dv.size() == 0) continue;
        double newval = *std::min_element(dv.begin(), dv.end());
        if (curval > newval) {
          curval = newval;
        }
      }
*/
    }
  }
  return curval;
}

double ClADataMap::getMaxValue() const {
  bool found = false;
  double curval = 0;
  for (const_iterator it = begin();it != end(); ++it) {
    if (it->second.type() == typeid(double)) {
      if (!found) {
        found = true;
        curval = anyc<double>(it->second);
      }
      else if (curval < anyc<double>(it->second)) {
        curval = anyc<double>(it->second);
      }
    }
    if (it->second.type() == typeid(ClCValue) && anyc<ClCValue>(it->second).hasValue<std::vector<double>>(VALUE)) {
      std::vector<double> dv = anyc<ClCValue>(it->second).value<std::vector<double>>(VALUE);
      if (dv.size() == 0) continue;
      // SJM 2023-01-06: exclude undefined values from check
      for (auto v : dv) {
        if (v != -999.25) {
          if (!found || v > curval) {
            curval = v;
            found = true;
          }
        }
      }
/*
      if (!found) {
        std::vector<double> dv = anyc<ClCValue>(it->second).value<std::vector<double>>(VALUE);
        if (dv.size() == 0)continue;
        found = true;
        if (dv.size() > 0) curval = *std::max_element(dv.begin(), dv.end());
      }
      else {
        std::vector<double> dv = anyc<ClCValue>(it->second).value<std::vector<double>>(VALUE);
        if (dv.size() == 0) continue;
        double newval = *std::max_element(dv.begin(), dv.end());
        if (curval < newval) {
          curval = newval;
        }
      }
*/
    }
  }
  return curval;
}


ClADataMap_PT ClADataMap::getBidirValuesInInterval(double top, double bottom) {
  ClADataMap_PT res(new ClADataMap(CDT_Bidir_Double_Point, EcUuid::sequential()));
  auto keys = getKeysInInterval(top, bottom, false);
  for (auto it=keys.begin(); it!=keys.end(); ++it) {
    if (*it==top) {
      auto v = getValueDownwards(top);
      if (!v.empty()) res->insert(top, v);
    } else if (*it==bottom) {
      auto v = getValueUpwards(bottom);
      if (!v.empty()) res->insert(bottom, v);

    } else {
      res->insertMultiList(*it, *it, values(*it));
    }
  }
  return res;
}

std::set<double> ClADataMap::getKeysInInterval(double topp, double bottomp, bool includefirstoutside) const {
  double top = capPrecision(topp);
  double bottom = capPrecision(bottomp);
  std::set<double> res;
  ClADataMap::const_iterator it = lowerBound(top);
  if (it != begin() && includefirstoutside) --it;
  for (; it != end() && it->first <= bottom; ++it) {
    res.insert(it->first);
  }
  if (it != end() && includefirstoutside) res.insert(it->first);
  return res;
}

ClADataMap::iterator ClADataMap::findMulti(const double& key, int ix) {
  if (_data.get<tag_depth>().count(capPrecision(key)) <= ix) return _data.get<tag_depth>().end();
  auto it = _data.get<tag_depth>().find(capPrecision(key));
  return boost::next(it, ix);
}

ClADataMap::iterator ClADataMap::replace(const double& top, const double &bottom, const boost::any& value) {
  iterator it = _data.get<tag_depth>().find(capPrecision(top));
  if (it == _data.get<tag_depth>().end()) {
    return _data.get<tag_depth>().insert(it, { top, bottom, value, EcUuid::sequential() });
  }
  else {
    replace(it, value);
    return it;
  }
}
ClADataMap::iterator ClADataMap::replace(const double& top, const boost::any& value) {
  return replace(top, top, value);
}

ClADataMap::iterator ClADataMap::find(const double& key) {
  return _data.get<tag_depth>().find(capPrecision(key));
}

ClADataMap::const_iterator ClADataMap::find(const double& key) const {
  return _data.get<tag_depth>().find(capPrecision(key));
}

ClADataMap::iterator ClADataMap::lowerBound(const double &key) {
  return _data.get<tag_depth>().lower_bound(capPrecision(key));
}

ClADataMap::const_iterator ClADataMap::lowerBound(const double &key) const {
  return _data.get<tag_depth>().lower_bound(key);
}

ClADataMap::iterator ClADataMap::upperBound(const double &key) {
  return _data.get<tag_depth>().upper_bound(capPrecision(key));
}

ClADataMap::const_iterator ClADataMap::upperBound(const double &key) const {
  return _data.get<tag_depth>().upper_bound(key);
}

QList<double> ClADataMap::keys() const {
  QList<double> res;
  for (auto it = begin(); it != end(); ++it) res.push_back(it->first);
  return res;
}

QList<double> ClADataMap::uniqueKeys() const {
  QSet<double> uniquekeys;
  for (auto it = begin(); it != end(); ++it) uniquekeys.insert(it->first);
  return uniquekeys.values();
}

QList<boost::any> ClADataMap::values() const {
  QList<boost::any> res;
  for (auto it = begin(); it != end(); ++it) res.push_back(it->second);
  return res;
}

QList<boost::any> ClADataMap::values(const double &key) const {
  QList<boost::any> res;
  for (auto it = _data.get<tag_depth>().find(key); it != end() && it->first == capPrecision(key); ++it) res.push_back(it->second);
  return res;
}

QList<DEntry> ClADataMap::entries(const double& key) const {
  QList<DEntry> res;
  for (auto it = _data.get<tag_depth>().find(key); it != end() && it->first == capPrecision(key); ++it) res.push_back(*it);
  return res;
}

std::set<DEntry> ClADataMap::entriesSet(const double& key) const {
  std::set<DEntry> res;
  for (auto it = _data.get<tag_depth>().find(key); it != end() && it->first == capPrecision(key); ++it) res.insert(*it);
  return res;
}

const boost::any ClADataMap::value(const double &key, const boost::any& def) const {
  auto it=_data.get<tag_depth>().find(capPrecision(key));
  if (it!=_data.get<tag_depth>().end()) {
    return it->second;
  }
  else {
    return def;
  }
}

std::vector<boost::any> ClADataMap::setSingleValue(const double &key, const boost::any & value) {
  std::vector<boost::any> res;
  for (auto it = find(key); it != end() && it->first == key; ++it) {
    res.push_back(it->second);
  }
  if (_data.get<tag_depth>().count(key) > 2) {
    EC_DBG_THROW();
  }
  _data.get<tag_depth>().erase(key);
  replace(key, value);
  return res;
}

std::vector<boost::any> ClADataMap::setValueUpwards(const double &key, const boost::any& value) {
  std::vector<boost::any> res;
  for (auto it = find(key); it != end() && it->first == key; ++it) {
    res.push_back(it->second);
  }
  if (_data.get<tag_depth>().count(key) == 0) {
    if (!value.empty()) {
      _data.get<tag_depth>().insert({ key, key, value, EcUuid::sequential() });
      _data.get<tag_depth>().insert({ key, key, boost::any(), EcUuid::sequential() });
    }
    return res;
  }
  if (_data.get<tag_depth>().count(key) == 1) {
    _data.get<tag_depth>().insert({ key, key, _data.get<tag_depth>().find(key)->second, EcUuid::sequential() });
    auto it = _data.get<tag_depth>().find(key);
    replace(it, value);
    return res;
  }
  if (_data.get<tag_depth>().count(key) == 2) {
    if (value.empty() && boost::next(_data.get<tag_depth>().find(key), 1)->second.empty()) {
      remove(key);
    }
    else {
      auto it = find(key);
      replace(it, value);
    }
    return res;
  }
  EC_DBG_THROW();
  return res;
}
/*
ClDataModifications ClADataMap::setValueUpwardsDM(const double &key, const boost::any& value, ClDocument_PT cldoc) {
  std::vector<ClDataModification> mod;
  for (auto it = find(key); it != end() && it->first == key; ++it) {
    res.push_back(it->second);
  }
  if (_data.count(key) == 0) {
    _data.insert(std::make_pair(key, value));
    _data.insert(std::make_pair(key, boost::any()));
//    return res;
  }
  if (_data.count(key) == 1) {
    _data.insert(std::make_pair(key, _data.find(key)->second));
    _data.find(key)->second = value;
    return res;
  }
  if (_data.count(key) == 2) {
    if (value.empty() && boost::next(_data.find(key), 1)->second.empty()) {
      remove(key);
    }
    else {
      find(key)->second = value;
    }
    return res;
  }
  EC_DBG_THROW();
}
*/
std::vector<boost::any> ClADataMap::setValueDownwards(const double &key, const boost::any& value) {
  std::vector<boost::any> res;
  for (auto it = find(key); it != end() && it->first == key; ++it) {
    res.push_back(it->second);
  }
  if (_data.get<tag_depth>().count(key) == 0) {
    if (!value.empty()) {
      _data.get<tag_depth>().insert({ key, key, boost::any(), EcUuid::sequential() });
      _data.get<tag_depth>().insert({ key, key, value, EcUuid::sequential() });
    }
    return res;
  }
  if (_data.get<tag_depth>().count(key) == 1) {
    _data.get<tag_depth>().insert({ key, key, value, EcUuid::sequential() });
    return res;
  }
  if (_data.get<tag_depth>().count(key) == 2) {
    if (value.empty() && _data.get<tag_depth>().find(key)->second.empty()) {
      remove(key);
    }
    else {
      auto it = boost::next(_data.get<tag_depth>().find(key), 1);
      replace(it, value);
    }
    return res;
  }
  EC_DBG_THROW();
  return res;
}
std::vector<boost::any> ClADataMap::setValueDir(const double &key, const boost::any& value, bool upwards) {
  if (upwards) {
    return setValueUpwards(key, value);
  }
  else {
    return setValueDownwards(key, value);
  }
}

boost::any ClADataMap::getValueUpwards(const double &key) const {
  auto it=find(key);
  if (it==end()) return boost::any();
  auto ti=boost::next(it,1);
  if (ti==end() || it->first != ti->first) return it->second;
  auto ti2 = boost::next(ti, 1);
  if (ti2 == end() || it->first != ti2->first) return it->second;
//  qDebug() << "getValueUpwards: More than 2 values in dual value data: Depth=" << QString::number(key/72/12, 'g', 18)+"ft/" +QString::number(key,'g',18) << " cnt=" << _data.count(key) << "V1=" << anyc<double>(values(key)[0]) << "V2=" << anyc<double>(values(key)[1]) << "V3=" << anyc<double>(values(key)[2]);
  EC_DBG_THROW();
  return it->second;
}

boost::any ClADataMap::getValueDownwards(const double &key) const {
  auto it = find(key);
  if (it==end()) return boost::any();
  auto ti = boost::next(it, 1);
  if (ti == end() || it->first != ti->first) return it->second;
  auto ti2 = boost::next(ti, 1);
  if (ti2 == end() || it->first != ti2->first) return ti->second;
//  qDebug() << "getValueDownwards: More than 2 values in dual value data: Depth=" << QString::number(key / 72 / 12, 'g', 18) + "ft/" + QString::number( key,'g', 18) << " cnt=" << _data.count(key) << "V1=" << anyc<double>(values(key)[0]) << "V2=" << anyc<double>(values(key)[1]) << "V3=" << anyc<double>(values(key)[2]);
  EC_DBG_THROW();
  return ti->second;
}
boost::any ClADataMap::getValueDir(const double &key, bool upwards)const {
  if (upwards) {
    return getValueUpwards(key);
  }
  else {
    return getValueDownwards(key);
  }
}


ClADataMap::iterator ClADataMap::erase(iterator it) { 
  return _data.get<tag_depth>().erase(it);
}

size_t ClADataMap::remove(const double &key) { 
  return _data.get<tag_depth>().erase(capPrecision(key));
}

void ClADataMap::insertAtEnd(const double& key, boost::any&& data) {
  _data.get<tag_depth>().insert(_data.get<tag_depth>().begin(), { capPrecision(key), capPrecision(key), std::move(data), EcUuid::sequential() });
}

ClADataMap::iterator ClADataMap::insert(const DEntry& entry) {
  return _data.get<tag_depth>().insert(entry).first;
}

ClADataMap::iterator ClADataMap::insert(const double& key, boost::any&& data) {
  return _data.get<tag_depth>().insert({ capPrecision(key), capPrecision(key), std::move(data), EcUuid::sequential() }).first;
}

ClADataMap::iterator ClADataMap::insert(const double &key, const boost::any& data) { 
  return _data.get<tag_depth>().insert({ capPrecision(key), capPrecision(key), data, EcUuid::sequential() }).first;
}


void ClADataMap::clear() { 
  _data.clear();
}

std::vector<boost::any> ClADataMap::replaceValues(const double &key, const std::vector<boost::any>& cals) {
  std::vector<boost::any> res;
  for (auto it = find(key); it != end() && it->first == key; ++it) {
    res.push_back(it->second);
  }
  remove(key);
  for (auto it = cals.begin(); it != cals.end(); ++it) {
    insertMulti(key, key, *it);
  }
  return res;
}

void ClADataMap::insertDEntryList(const QList<DEntry>& data) {
  if (data.empty()) return;
  auto dit = data.rbegin();
  ClADataMap::iterator it = _data.get<tag_depth>().insert(*dit).first;

  for (; dit != data.rend(); ++dit) {
    DEntry de = *dit;
    it = _data.get<tag_depth>().insert(it,de);
  }
}

void ClADataMap::insertMultiList(const double &top, const double& bot, const QList<boost::any>& data) {
  for (auto it = data.begin(); it != data.end();++it) {
    _data.insert({ capPrecision(top), capPrecision(bot), *it , EcUuid::sequential() });
  }
}

size_t ClADataMap::insertMulti(const double &top, const double &bot, const boost::any& data, int ix) {
  if ((int)_data.get<tag_depth>().count(top) < ix) ix = (int)_data.get<tag_depth>().count(top);

  _data.get<tag_depth>().insert(boost::next(_data.get<tag_depth>().find(top), ix), { capPrecision(top), capPrecision(bot), data, EcUuid::sequential() });
  return ix;
}

size_t ClADataMap::insertMulti(const double &top, const double& bot, const boost::any& data) {
  _data.get<tag_depth>().insert({ capPrecision(top), capPrecision(bot), data, EcUuid::sequential() });
  return _data.get<tag_depth>().count(capPrecision(top)) - 1;
}

double ClADataMap::getDouble(double depthp) const {
  double depth = capPrecision(depthp);
  if (!contains(depth) || find(depth)->second.type() != typeid(double)) {
#ifdef DEBUG
    qDebug() << "getDouble used on missing or incorrect value";
#endif
    return 0;
  }
  return anyc<double>(find(depth)->second);
}

int ClADataMap::getInt(double depthp) const {
  double depth = capPrecision(depthp);
  if (!contains(depth) || find(depth)->second.type() != typeid(int)) {
#ifdef DEBUG
    qDebug() << "getInt used on missing or incorrect value";
#endif
    return 0;
  }
  return anyc<int>(find(depth)->second);
}

ClAGenericMap_PT ClADataMap::getCopyNewIds() const {
  ClADataMap_PT res(new ClADataMap(_dataType, EcUuid::sequential()));
  res->_unitName = _unitName;
  for (auto it = begin(); it != end(); ++it) {
    DEntry ent= *it;
    ent.uuid = EcUuid::sequential();
    res->insert(ent);
  }
  return res;
}


QColor ClADataMap::getColor(double depthp) const {
  double depth = capPrecision(depthp);
  if (!contains(depth) || find(depth)->second.type() != typeid(QColor)) {
#ifdef DEBUG
    qDebug() << "getColor used on missing or incorrect value";
#endif
    return QColor();
  }
  return anyc<QColor>(find(depth)->second);
}

double ClADataMap::getTopDepth() const {
  if (empty()) return 0;
  return begin()->first;
}

double ClADataMap::getBottomDepth() const {
  if (empty()) return 0;
  return _data.get<tag_bottom>().rbegin()->bottom;

}
ClADataMap::iterator ClADataMap::getFirstIteratorAbove(double vp) {
  double v = capPrecision(vp);
  iterator it = upperBound(v);
  if (it == begin()) return it;
  return boost::prior(it, 1);
}

ClADataMap::const_iterator ClADataMap::getFirstIteratorAbove(double vp) const {
  double v = capPrecision(vp);
  const_iterator it = upperBound(v);
  if (it == begin()) return it;
  return boost::prior(it, 1);
}

double ClADataMap::getDualDataPosAbove(double vp, bool *valid, bool exclusive) const {
  if (!valid) {
    qDebug() << "getDualDataPosAbove must have valid parameter";
  }
  double above = getDataPosAbove(vp, valid, exclusive);
  if (!*valid) return above;
  auto it = find(above);
  for (; count(it->first) <= 1; --it) {
    if (it == begin()) {
      *valid = false;
      return above;
    }
  }
  return it->first;
}

double ClADataMap::getDualDataPosBelow(double vp, bool *valid, bool exclusive) const {
  if (!valid) {
    qDebug() << "getDualDataPosBelow must have valid parameter";
  }
  double below = getDataPosBelow(vp, valid, exclusive);
  if (!*valid) return below;
  auto it = find(below);
  for (; count(it->first) <= 1; ++it) {
    if (boost::next(it,1) == end()) {
      *valid = false;
      return below;
    }
  }

  return it->first;
}

double ClADataMap::getDataPosAbove(double vp, bool *valid, bool exclusive) const {
  double v = capPrecision(vp);

  if (valid) *valid = false;
  if (!count()) return 0;
  if (exclusive) {
    auto it = find(v);
    if (it != end()) {
      if (v == begin()->first) {
        return 0;
      }
      else {
        if (valid) *valid = true;
        return boost::prior(it, 1)->first;
      }
    }
  }

  const_iterator it = upperBound(v);
  if (it == begin()) {
    return 0;
  }
  --it;
  if (valid) *valid = true;
  return it->first;
}

double ClADataMap::getDataPosAboveBottom(double vp, bool* valid, bool exclusive) const {
  double v = capPrecision(vp);

  if (valid) *valid = false;
  if (!count()) return 0;
  if (exclusive) {
    auto it = find(v);
    if (it != end()) {
      if (v == begin()->first) {
        return 0;
      }
      else {
        if (valid) *valid = true;
        return boost::prior(it, 1)->bottom;
      }
    }
  }

  const_iterator it = upperBound(v);
  if (it == begin()) {
    return 0;
  }
  --it;
  if (valid) *valid = true;
  return it->bottom;
}

// SJM 2021-02-23: added function to just get bottom value of entry at depth
double ClADataMap::getDataPosBottom(double vp, bool* valid) const {
  double v = capPrecision(vp);

  if (valid) *valid = false;
  if (!count()) return 0;
  auto it = find(v);
  if (it != end()) {
    if (valid) *valid = true;
    return it->bottom;
  }

  return 0;
}

double ClADataMap::getDataPosBelow(double vp, bool *valid, bool exclusive) const {
  double v = capPrecision(vp);
  if (valid) *valid = false;
  if (!count()) return 0;
  if (exclusive) {
    auto it = find(v);
    if (it != end()) {
      if (v == boost::prior(end(), 1)->first) {
        return 0;
      }
      else {
        if (valid) *valid = true;
        while (it->first == v) ++it;
        return it->first;
      }
    }
  }
  const_iterator it = lowerBound(v);
  if (it != end()) {
    if (valid) *valid = true;
    return it->first;
  }
  else {
    return 0;
  }
}

double ClADataMap::getClosestIntervalBorder(double vp, bool* ok) const {
//  double v = capPrecision(vp);
  if (empty()) {
    if (ok) *ok = false;
    return 0;
  }
  if (ok) *ok = true;
  std::vector<double> candidates;
  double bestcandidate = 0;
  double bestcandidatedist = -1;
  auto it = _data.get<tag_depth>().lower_bound(vp);
  if (it != _data.get<tag_depth>().end()) {
    bestcandidate = it->first;
    bestcandidatedist = fabs(it->first - vp);
  }
  if (it != _data.get<tag_depth>().begin()) {
    --it;
    double newdist = fabs(it->first - vp);
    if (newdist < bestcandidatedist || bestcandidatedist < 0) {
      bestcandidatedist = newdist;
      bestcandidate = it->first;
    }
  }
  auto ti = _data.get<tag_bottom>().lower_bound(vp);
  if (ti != _data.get<tag_bottom>().end()) {
    double newdist = fabs(ti->bottom - vp);
    if (newdist < bestcandidatedist || bestcandidatedist<0) {
      bestcandidatedist = newdist;
      bestcandidate = ti->bottom;
    }
  }
  if (ti != _data.get<tag_bottom>().begin()) {
    --ti;
    double newdist = fabs(ti->bottom - vp);
    if (newdist < bestcandidatedist || bestcandidatedist < 0) {
      bestcandidatedist = newdist;
      bestcandidate = it->bottom;
    }
  }
  return bestcandidate;
}

double ClADataMap::getClosestDataPos(double vp, bool *ok) const {
  double v = capPrecision(vp);
  if (ok) *ok=false;
  if (!count()) return 0;
  if (ok) *ok=true;
  const_iterator it = lowerBound(v);
  if (it == end()) {
    --it;
    return it->first;
  }
  if (it == begin()) {
    return it->first;
  }
  double v1 = it->first;
  --it;
  double v2 = it->first;
  if (fabs(v1 - v) < fabs(v2 - v)) {
    return v1;
  }
  else {
    return v2;
  }
}

void ClADataMap::rwXml(const QString &key, QDomElement &dom, bool isRead, bool isArchive, bool ignoreimages, const QMap<QString, QString> &imagemapnames) {
  EC_UNUSED(imagemapnames, ignoreimages)
  (void)isArchive;
  QDomElement elem;
  if (isRead) {
    elem = dom.firstChildElement(key);
  }
  else {
    elem = dom.ownerDocument().createElement("Map");
    elem.setAttribute("name", key);
    QBuffer buf;
    QTextStream ts(&buf);
    for (ClADataMap::iterator it = begin(); it != end(); ++it) {
      ts << it->first << " ";
      ts << boostAnyToString(it->second) << "\n";
    }
    dom.appendChild(elem);
  }
}

double ClADataMap::getBottomExtent() const {
  if (empty()) return 0;
  if (dataType() == CDT_Notes) {
    return boost::prior(end(), 1)->bottom;
  } else if (dataType() == CDT_Image) {
    return boost::prior(end(), 1)->bottom;
  }
  else {
    return boost::prior(end(), 1)->first;
  }
}

double ClADataMap::getTopExtent() const {
  if (empty()) return 0;
  return begin()->first;
}

void ClADataMap::purgeImageData() {
  for (auto it = begin(); it != end(); ++it) {
    if (it->second.type() == typeid(ClColumnImage)) {
      const ClColumnImage &clim = boost::any_cast<const ClColumnImage&>(it->second);
      if (clim._image->sourceType!=Ec::SourceDatabase) clim._image->unloadMetafile();
    }
  }
}


void ClADataMap::reloadImageData() {
  for (auto it = begin(); it != end(); ++it) {
    if (it->second.type() == typeid(ClColumnImage)) {
      const ClColumnImage &clim = boost::any_cast<const ClColumnImage&>(it->second);
      if (clim._image->sourceType != Ec::SourceDatabase) {
        clim._image->loadImage(true);
        clim._image->generateRasterMetafile();
      }
    }
  }
}

double ClADataMap::getDoubleInterpolated(double depth, bool *ok, bool preferdownward) const {

  if (ok) *ok = true;
  auto it = find(depth);
  if (it != end()) {
    auto nxtit = boost::next(it, 1);
    if (nxtit != end() && nxtit->first == it->first) {
      if (preferdownward) {
        if (!nxtit->second.empty()) {
          return anyc<double>(nxtit->second);
        }
        else {
          if (ok) *ok = false;
          return 0;
        }
      }
      else {
        if (!it->second.empty()) {
          return anyc<double>(it->second);
        }
        else {
          if (ok) *ok = false;
          return 0;
        }
      }
    }
    if (!it->second.empty())  return anyc<double>(it->second);
  }
  bool valid;
  double top = getDataPosAbove(depth, &valid, true);
  if (!valid) {
    if (ok) *ok = false;
    return 0;
  }
  boost::any atv = getValueDownwards(top);
  double tv = 0;
  if (!atv.empty()) tv=anyc<double>(atv);
//  if (contains(top)) return tv;
  double bottom = getDataPosBelow(depth, &valid, true);
  if (!valid) {
    if (ok) *ok = false;
    return 0;
  }
  boost::any abv = getValueUpwards(bottom);
  double bv = 0;
  if (!abv.empty()) bv = anyc<double>(abv);
//  if (contains(bottom)) return bv;
  if (abv.empty() || atv.empty()) {
    if (ok) *ok = false;
    return 0;
  }
  return tv + (bv - tv)*(depth - top) / (bottom - top);
}

std::vector<std::pair<double, double> > ClADataMap::getDataGapsInInterval(double top, double bot) const {
  if (top >= bot) return {};
  std::vector<std::pair<double, double> > res;
  ClADataMap intervals;
  ClDataMapRange rng = getRange(top, bot);
  if (rng.begin() == rng.end()) {
    res.push_back({ top,bot });
  }
  else {
    for (auto it = rng.begin(); it != rng.end(); ++it) {
      intervals.insert(*it);
    }
    double gapstart = top;
    for (auto it = intervals.begin(); it != intervals.end(); ++it) {
      if (gapstart < it->first) {
        res.push_back({ gapstart, MIN(it->first, bot) });
      }
      gapstart = MAX(gapstart, it->bottom);
    }
  }
  return res;
}


double ClADataMap::getFirstEmptyDepthBelow(double dpth) const {
  double bot = getBottomDepth();
  if (dpth > bot) return dpth;
  auto ivs = getDataGapsInInterval(dpth, bot);
  if (ivs.empty()) return bot;
  return ivs.front().first;
}

double ClADataMap::getFirstEmptyDepthAbove(double dpth) const {
  double top = getTopDepth();
  if (dpth < top) return dpth;
  auto ivs = getDataGapsInInterval(top, dpth);
  if (ivs.empty()) return top;
  return ivs.back().second;
}

ClAGenericMap_PT ClAGenericMap::createMap(ClDataType dt) {
  if (dt == CDT_Depth_Map) {
    return ClDepthCorrection_PT(new ClDepthCorrection());
  }
  else {
    return ClADataMap_PT(new ClADataMap());
  }
}

ClADataMap_PT ClADataMap::extractInterval(double top, double bottom) const {
  ClADataMap_PT res(new ClADataMap(_dataType, EcUuid::sequential()));
  if (_dataType == CDT_Bidir_Double_Point) {
    if (count(top) >= 1) res->insert(top, getValueDownwards(top));
    if (count(bottom) >= 1) res->insert(bottom, getValueUpwards(bottom));
    for (auto it = getFirstIteratorAbove(top); it != end() && it->first < bottom; ++it) {
      if (it->first <= top) continue;
      res->insertMultiList(it->first, it->first, values(it->first));
      while (it != end() && boost::next(it,1)!=end() && boost::next(it, 1)->first == it->first) ++it;
    }
  }
  return res;
}

/*
bool ClADataMap::getTightDataRange(double top, double bottom, double* maxtop, double* maxbottom) {
  if (_dataType == CDT_Int_Interval || _dataType == CDT_Color_Interval || _dataType == CDT_Age_Interval || _dataType == CDT_Lithology_Interval || _dataType == CDT_String_Interval) {
    bool avalidtop, bvalidtop, avalidbot;
    double abovetop = getDataPosAbove(top, &avalidtop, false);
    double belowtop = getDataPosBelow(top, &bvalidtop, true);
    double abovebot = getDataPosAbove(top, &avalidbot, true);
    if (!avalidtop && !bvalidtop) return false;
    if (avalidtop && !bvalidtop) return false;
    if (maxtop) {
      if (!avalidtop && bvalidtop) *maxtop = belowtop;
      if (avalidtop && bvalidtop) {
        if (find(abovetop)->second.empty()) {
          *maxtop = belowtop;
        }
        else {
          *maxtop = top;
        }
      }
    }
    if (maxbottom) {
      if (avalidbot && find(abovebot)->second.empty()) {
        *maxbottom = abovebot;

      }
      else {
        *maxbottom = bottom;
      }
    }
  }
  else if (_dataType == CDT_Notes) {
    bool found = false;
    for (auto it = begin(); it != end(); ++it) {
      ClNotesEntry ne = anyc<ClNotesEntry>(it->second);
      if (it->first<bottom && ne.bottomdepth>top) {
//        if (maxtop) *maxtop=
      }
//      if (it->second.)
    }
  }
}*/

//ClIntervalMap *_ivmap;
//ClADataMap::iterator dmit;

class ClCValue_equals : public boost::static_visitor<bool> {
public:
  template <typename T, typename U>
  bool operator()(const T&, const U&) const
  {
    return false; // cannot compare different types
  }

  template <typename T>
  bool operator()(const T& lhs, const T& rhs) const
  {
    return lhs == rhs;
  }

  template<class T>
  boost::optional<bool> checkTypeAny(const boost::any& a1, const boost::any& a2) const {
    if (a1.type() != typeid(T)) return {};
    return boost::any_cast<T>(a1) == boost::any_cast<T>(a2);
  }
  bool operator()(const QMap<QString, boost::any>& lhs, const QMap<QString, boost::any>& rhs) const
  {
    const QMap<QString, boost::any>* rhp = boost::get<const QMap<QString, boost::any>*>(&rhs);
    const QMap<QString, boost::any>* lhp = boost::get<const QMap<QString, boost::any>*>(&lhs);
    if (lhp->size() != rhp->size()) return false;
    for (auto it = lhp->begin(); it != lhp->end(); ++it) {
      auto ti = rhp->find(it.key());
      if (ti == rhp->end()) return false;
      if (ti->empty()!= it->empty()) return false;
      if (ti->type() != it->type()) return false;
      boost::optional<bool> res;
      res = checkTypeAny<double>(it.value(), ti.value());
      if (res && !*res) return false;
      res = checkTypeAny<int>(it.value(), ti.value());
      if (res && !*res) return false;
      res = checkTypeAny<QPolygonF>(it.value(), ti.value());
      if (res && !*res) return false;
      res = checkTypeAny<ClDepthPointVector>(it.value(), ti.value());
      if (res && !*res) return false;
      res = checkTypeAny<QString>(it.value(), ti.value());
      if (res && !*res) return false;
      res = checkTypeAny<ClColumnImage>(it.value(), ti.value());
      if (res && !*res) return false;
      res = checkTypeAny<QColor>(it.value(), ti.value());
      if (res && !*res) return false;
    }
    return true;
  }
};

ClCValue::ClCValue() {
}

ClCValue::ClCValue(int v) {
  _vmap = v;
}
ClCValue::ClCValue(double v) {
  _vmap = v;
}
ClCValue::ClCValue(QString v) {
  _vmap = v;
}

bool ClCValue::operator==(const ClCValue& other) const {
  return boost::apply_visitor(ClCValue_equals(), this->_vmap, other._vmap);
}

bool  ClCValue::hasValueAny(const QString& str) const { 
  switch (_vmap.which()) {
    case 0: return false; //boost::blank
    case 1: return boost::get<QMap<QString, boost::any>>(&_vmap)->contains(str);
    default:
      return str == "VALUE";
      break;
  }
}

bool ClCValue::hasPointSet() const {
  if (hasValue<ClDepthPointVector>(POINTSET_REL)) return true;
  if (hasValue<ClDepthPointVector>(POINTSET_ABS)) return true;
  if (hasValue<ClDepthPointVector>(POINTSET_SCL)) return true;
  if (hasValue<ClDepthPointVector>(POINTSET_PHI)) return true;
  if (hasValue<ClDepthPointVector>(POINTSET_XPIC)) return true;
  return false;
}

bool ClCValue::hasPosition() const {
  if (hasValue<double>(POSITION_REL)) return true;
  if (hasValue<double>(POSITION_ABS)) return true;
  if (hasValue<double>(POSITION_SCL)) return true;
  if (hasValue<double>(POSITION_PHI)) return true;
  if (hasValue<double>(POSITION_XPIC)) return true;
  if (hasValue<double>(POSITION_YPIC)) return true;
  return false;

}

double ClCValue::pointSetTopDepth() const {
  if (hasValue<ClDepthPointVector>(POINTSET_PHI)) return value<ClDepthPointVector>(POINTSET_PHI).topDepth();
  if (hasValue<ClDepthPointVector>(POINTSET_SCL)) return value<ClDepthPointVector>(POINTSET_SCL).topDepth();
  if (hasValue<ClDepthPointVector>(POINTSET_REL)) return value<ClDepthPointVector>(POINTSET_REL).topDepth();
  if (hasValue<ClDepthPointVector>(POINTSET_ABS)) return value<ClDepthPointVector>(POINTSET_ABS).topDepth();
  return 0;
}

double ClCValue::pointSetBottomDepth() const {
  if (hasValue<ClDepthPointVector>(POINTSET_PHI)) return value<ClDepthPointVector>(POINTSET_PHI).bottomDepth();
  if (hasValue<ClDepthPointVector>(POINTSET_SCL)) return value<ClDepthPointVector>(POINTSET_SCL).bottomDepth();
  if (hasValue<ClDepthPointVector>(POINTSET_REL)) return value<ClDepthPointVector>(POINTSET_REL).bottomDepth();
  if (hasValue<ClDepthPointVector>(POINTSET_ABS)) return value<ClDepthPointVector>(POINTSET_ABS).bottomDepth();
  return 0;
}

template<class T> void ClCValue::dbgCheckType(const QString& str) const{
  // AVA: sanity check in debug version to ensure that we have a type/string match on builtin properties.
  if ((str == POSITION_ABS || str == POSITION_PHI || str == POSITION_REL || str == POSITION_SCL || str == POSITION_XPIC || str == POSITION_YPIC) && strcmp(typeid(T).name(), "double") != 0) EC_THROW("Bad check on builtin in ClCValue::hasValue", str, typeid(T).name());
  if ((str == POINTSET_ABS || str == POINTSET_PHI || str == POINTSET_REL || str == POINTSET_SCL) && strcmp(typeid(T).name(), "class ClDepthPointVector") != 0) EC_THROW("Bad check on builtin in ClCValue::hasValue", str, typeid(T).name());
}

template<class T> 
bool  ClCValue::hasValue(const QString& str) const { 
#ifdef DEBUG
  dbgCheckType<T>(str);
#endif
  if (_vmap.which() == 0) return false;
  if (str == "VALUE") {
    if (_vmap.which() == 2) return typeid(T)==typeid(int);
    if (_vmap.which() == 3) return typeid(T)==typeid(double);
    if (_vmap.which() == 4) return typeid(T)==typeid(QString);
    if (_vmap.which() == 1) {
      return boost::get<QMap<QString, boost::any> >(&_vmap)->value(str).type() == typeid(T);
    }
    return false;
  }
  else {
    if (_vmap.which() != 1) return false;
    return boost::get<QMap<QString, boost::any> >(&_vmap)->value(str).type() == typeid(T);
  }
}

template bool ClCValue::hasValue<QPolygonF>(const QString& str) const;
template bool ClCValue::hasValue<bool>(const QString& str) const;
template bool ClCValue::hasValue<double>(const QString& str) const;
template bool ClCValue::hasValue<int>(const QString& str) const;
template bool ClCValue::hasValue<QString>(const QString& str) const;
template bool ClCValue::hasValue<ClColumnImage>(const QString& str) const;
//template bool ClCValue::hasValue<ClRefDepthPosition>(const QString& str) const;
//template bool ClCValue::hasValue<ClRefDepthPointVector>(const QString& str) const;
template bool ClCValue::hasValue<ClDepthPointVector>(const QString& str) const;
template bool ClCValue::hasValue<std::vector<double>>(const QString& str) const;

boost::any  ClCValue::valueAny(const QString& str) const { 
  switch (_vmap.which()) {
    default:
    case 0: return {}; //boost::blank
    case 1: return boost::get<QMap<QString, boost::any>>(&_vmap)->value(str);;
    case 2: if (str == "VALUE") return boost::any(boost::get<int>(_vmap)); else return {};
    case 3: if (str == "VALUE") return boost::any(boost::get<double>(_vmap)); else return {};
    case 4: if (str == "VALUE") return boost::any(boost::get<QString>(_vmap)); else return {};
  }
  return {};
}


template<class T> 
T ClCValue::value(const QString& str) const { 
#ifdef DEBUG
  dbgCheckType<T>(str);
#endif
  switch (_vmap.which()) {
    default:
    case 0: return {}; //boost::blank
    case 1: {
              auto mp = boost::get<QMap<QString, boost::any>>(&_vmap);
              auto it = mp->find(str);
              if (it == mp->end() || it->type()!=typeid(T)) return {};
              return boost::any_cast<T>(*it);
            }
  }
  return {};
}

template<>
int ClCValue::value(const QString& str) const {
#ifdef DEBUG
  dbgCheckType<int>(str);
#endif
  switch (_vmap.which()) {
  default:
  case 0: return {}; //boost::blank
  case 1: {
    auto mp = boost::get<QMap<QString, boost::any>>(&_vmap);
    auto it = mp->find(str);
    if (it == mp->end() || it->type() != typeid(int)) return {};
    return boost::any_cast<int>(*it);
  }
  case 2: if (str == "VALUE") return boost::get<int>(_vmap); else return {};
  }
  return {};
}
template<>
double ClCValue::value(const QString& str) const {
#ifdef DEBUG
  dbgCheckType<double>(str);
#endif
  switch (_vmap.which()) {
  default:
  case 0: return {}; //boost::blank
  case 1: {
    auto mp = boost::get<QMap<QString, boost::any>>(&_vmap);
    auto it = mp->find(str);
    if (it == mp->end() || it->type() != typeid(double)) return {};
    return boost::any_cast<double>(*it);
  }
  case 3: if (str == "VALUE") return boost::get<double>(_vmap); else return {};
  }
  return {};
}

template<>
QString ClCValue::value(const QString& str) const {
#ifdef DEBUG
  dbgCheckType<QString>(str);
#endif
  switch (_vmap.which()) {
  default:
  case 0: return {}; //boost::blank
  case 1: {
    auto mp = boost::get<QMap<QString, boost::any>>(&_vmap);
    auto it = mp->find(str);
    if (it == mp->end() || it->type() != typeid(QString)) return {};
    return boost::any_cast<QString>(*it);
  }
  case 4: if (str == "VALUE") return boost::get<QString>(_vmap); else return {};
  }
  return {};
}

template bool ClCValue::value<bool>(const QString& str) const;
template QPolygonF ClCValue::value<QPolygonF>(const QString& str) const;
//template ClRefDepthPosition ClCValue::value<ClRefDepthPosition>(const QString& str) const;
//template ClRefDepthPointVector ClCValue::value<ClRefDepthPointVector>(const QString& str) const;
template ClDepthPointVector ClCValue::value<ClDepthPointVector>(const QString& str) const;
template ClColumnImage ClCValue::value<ClColumnImage>(const QString& str) const;
template int ClCValue::value<int>(const QString& str) const;
template std::vector<double> ClCValue::value<std::vector<double>>(const QString& str) const;

template<class T> bool ClCValue::hasValueBase() const {
  if (_vmap.which() != 1) return false;
  auto mp = boost::get<const QMap<QString, boost::any>&>(_vmap);
  auto it = mp.find("VALUE");
  if (it == mp.end()) return false;
  return it->type() == typeid(T);
}
template<> bool ClCValue::hasValueBase<int>() const {
  if (_vmap.which() == 2) return true;
  if (_vmap.which() != 1) return false;
  auto mp = boost::get<const QMap<QString, boost::any>&>(_vmap);
  auto it = mp.find("VALUE");
  if (it == mp.end()) return false;
  return it->type() == typeid(int);
}
template<> bool ClCValue::hasValueBase<double>() const {
  if (_vmap.which() == 3) return true;
  if (_vmap.which() != 1) return false;
  auto mp = boost::get<const QMap<QString, boost::any>&>(_vmap);
  auto it = mp.find("VALUE");
  if (it == mp.end()) return false;
  return it->type() == typeid(double);
}
template<> bool ClCValue::hasValueBase<QString>() const {
  if (_vmap.which() == 4) return true;
  if (_vmap.which() != 1) return false;
  auto mp = boost::get<const QMap<QString, boost::any>&>(_vmap);
  auto it = mp.find("VALUE");
  if (it == mp.end()) return false;
  return it->type() == typeid(QString);
}

template bool ClCValue::hasValueBase<bool>() const;
template bool ClCValue::hasValueBase<QPolygonF>() const;
template bool ClCValue::hasValueBase<ClDepthPointVector>() const;
template bool ClCValue::hasValueBase<ClColumnImage>() const;
template bool ClCValue::hasValueBase<int>() const;
template bool ClCValue::hasValueBase<std::vector<double> >() const;

template<class T>
T ClCValue::valueBase() const {
  switch (_vmap.which()) {
  default:
    case 0: return {}; //boost::blank
    case 1: {
      auto mp = boost::get<QMap<QString, boost::any>>(&_vmap);
      auto it = mp->find("VALUE");
      if (it == mp->end() || it->type() != typeid(T)) return {};
      return boost::any_cast<T>(*it);
    }
  }
  return {};
}
template<>
int ClCValue::valueBase() const {
  switch (_vmap.which()) {
  default:
  case 0: return {}; //boost::blank
  case 1: {
    auto mp = boost::get<QMap<QString, boost::any>>(&_vmap);
    auto it = mp->find("VALUE");
    if (it == mp->end() || it->type() != typeid(int)) return {};
    return boost::any_cast<int>(*it);
  }
  case 2: return boost::get<int>(_vmap);
  }
  return {};
}
template<>
double ClCValue::valueBase() const {
  switch (_vmap.which()) {
  default:
  case 0: return {}; //boost::blank
  case 1: {
    auto mp = boost::get<QMap<QString, boost::any>>(&_vmap);
    auto it = mp->find("VALUE");
    if (it == mp->end() || it->type() != typeid(double)) return {};
    return boost::any_cast<double>(*it);
  }
  case 3: return boost::get<double>(_vmap);
  }
  return {};
}

template<>
QString ClCValue::valueBase() const {
  switch (_vmap.which()) {
  default:
  case 0: return {}; //boost::blank
  case 1: {
    auto mp = boost::get<QMap<QString, boost::any>>(&_vmap);
    auto it = mp->find("VALUE");
    if (it == mp->end() || it->type() != typeid(QString)) return {};
    return boost::any_cast<QString>(*it);
  }
  case 4: return boost::get<QString>(_vmap); 
  }
  return {};
}

template bool ClCValue::valueBase<bool>() const;
template QPolygonF ClCValue::valueBase<QPolygonF>() const;
template ClDepthPointVector ClCValue::valueBase<ClDepthPointVector>() const;
template ClColumnImage ClCValue::valueBase<ClColumnImage>() const;
template int ClCValue::valueBase<int>() const;
template std::vector<double> ClCValue::valueBase<std::vector<double> >() const;

void  ClCValue::setValueAny(const QString& str, const boost::any& v) { 
  bool isint = (v.type() == typeid(int));
  bool isdouble = (v.type() == typeid(double));
  bool isQString = (v.type() == typeid(QString));
  switch (_vmap.which()) {
    case 1: 
      if (boost::get<QMap<QString, boost::any>>(&_vmap)->size() == 1 && boost::get<QMap<QString, boost::any>>(&_vmap)->contains("VALUE") && str == "VALUE" && (isint || isdouble || isQString)) {
        if (isint) _vmap = boost::any_cast<int>(v);
        if (isdouble) _vmap = boost::any_cast<double>(v);
        if (isQString) _vmap = boost::any_cast<QString>(v);
      }
      else {
        boost::get<QMap<QString, boost::any>>(&_vmap)->insert(str, v);
      }
      break;
    default: {
      if (str == "VALUE" && (isint || isdouble || isQString)) {
        if (isint) _vmap = boost::any_cast<int>(v); 
        if (isdouble) _vmap = boost::any_cast<double>(v);
        if (isQString) _vmap = boost::any_cast<QString>(v);        
      }
      else {
        QMap<QString, boost::any> mp;
        mp.insert(str, v);
        if (_vmap.which()==2) mp.insert("VALUE", boost::get<int>(_vmap));
        if (_vmap.which()==3) mp.insert("VALUE", boost::get<double>(_vmap));
        if (_vmap.which()==4) mp.insert("VALUE", boost::get<QString>(_vmap));
        _vmap = mp;
      }
    }
  }
}

template<class T> void clCValueSetIfcompatible(const T&, ClCValue&) {};
template<> void clCValueSetIfcompatible<int>(const int& i, ClCValue& me) { me._vmap = i; };
template<> void clCValueSetIfcompatible<double>(const double& i, ClCValue& me) { me._vmap = i; };
template<> void clCValueSetIfcompatible<QString>(const QString& i, ClCValue& me) { me._vmap = i; };

template<class T> void  ClCValue::setValue(const QString& str, const T& v) { 
#ifdef DEBUG
  dbgCheckType<T>(str);
#endif
  bool isint = (typeid(T) == typeid(int));
  bool isdouble = (typeid(T) == typeid(double));
  bool isQString = (typeid(T) == typeid(QString));
  switch (_vmap.which()) {
  case 1:
    if (boost::get<QMap<QString, boost::any>>(&_vmap)->size() == 1 && boost::get<QMap<QString, boost::any>>(&_vmap)->contains("VALUE") && str == "VALUE" && (isint || isdouble || isQString)) {
      clCValueSetIfcompatible(v,*this);
    }
    else {
      boost::get<QMap<QString, boost::any>>(&_vmap)->insert(str, v);
    }
    break;
  default: {
      if (str == "VALUE" && (isint || isdouble || isQString)) {
        clCValueSetIfcompatible(v, *this);
      }
      else {
        QMap<QString, boost::any> mp;
        mp.insert(str, v);
        if (_vmap.which() == 2) mp.insert("VALUE", boost::get<int>(_vmap));
        if (_vmap.which() == 3) mp.insert("VALUE", boost::get<double>(_vmap));
        if (_vmap.which() == 4) mp.insert("VALUE", boost::get<QString>(_vmap));
        _vmap = mp;
      }
    }
  }
}

//template void ClCValue::setValue<ClRefDepthPosition>(const QString& str, const ClRefDepthPosition &v);
//template void ClCValue::setValue<ClRefDepthPointVector>(const QString& str, const ClRefDepthPointVector &v);
template void ClCValue::setValue<ClDepthPointVector>(const QString& str, const ClDepthPointVector &v);
template void ClCValue::setValue<QColor>(const QString& str, const QColor &v);
template void ClCValue::setValue<QPolygonF>(const QString& str, const QPolygonF &v);
template void ClCValue::setValue<QString>(const QString& str, const QString &v);
template void ClCValue::setValue<double>(const QString& str, const double &v);
template void ClCValue::setValue<int>(const QString& str, const int &v);
template void ClCValue::setValue<std::vector<double> >(const QString& str, const std::vector<double>&v);

ClCValue& ClCValue::operator=(int v) {
  _vmap = v;
  return *this;
}

ClCValue& ClCValue::operator=(double v) {
  _vmap = v;
  return *this;
}

ClCValue& ClCValue::operator=(const QString& v) {
  _vmap = v;
  return *this;
}

int ClCValue::count() const {
  switch (_vmap.which()) {
    case 0: return 0;
    case 1: return boost::get<const QMap<QString, boost::any>&>(_vmap).size();
    default: return 1;
  }
}

template<class T> bool ClCValue::isBaseValue() const {
  if (_vmap.which() != 1) return false;
  auto mp = boost::get<const QMap<QString, boost::any>&>(_vmap);
  if (mp.size() == 1 && mp.begin().key() == "VALUE" && mp.begin().value().type() == typeid(T)) return true;
  return false;
}

template<> bool ClCValue::isBaseValue<int>() const {
  return _vmap.which() == 2;
}

template<> bool ClCValue::isBaseValue<double>() const {
  return _vmap.which() == 3;
}

template<> bool ClCValue::isBaseValue<QString>() const {
  return _vmap.which() == 4;
}


boost::any ClCValue::getIndexedAny(int i) {
  switch (_vmap.which()) {
    default:
    case 0: return {};
    case 1: {
      auto mp = boost::get<const QMap<QString, boost::any>&>(_vmap);
      if (i < 0 || i >= mp.size()) return {};
      return boost::next(mp.begin(), i).value();
    }
    case 2: if (i == 0) return boost::get<int>(_vmap); else return {};
    case 3: if (i==0) return boost::get<double>(_vmap); else return {};
    case 4: if (i==0) return boost::get<QString>(_vmap); else return {};
  }
}


void  ClCValue::clearValue(const QString& str) { 
  switch (_vmap.which()) {
    case 0: return; //boost::blank
    case 1: boost::get<QMap<QString, boost::any>>(&_vmap)->remove(str); return;
    default:
      if (str == "VALUE") _vmap=boost::blank();
      break;
  }
}
bool  ClCValue::operator !=(const ClCValue& other) const {
  return !(*this == other); 
}

QMap<QString, boost::any> ClCValue::asMap() const {
  switch (_vmap.which()) {
    default:
    case 0: return {}; //boost::blank
    case 1: return boost::get<QMap<QString, boost::any> >(_vmap);
    case 2: return { {"VALUE",boost::any(boost::get<int>(_vmap))} };
    case 3: return { {"VALUE", boost::any(boost::get<double>(_vmap))} }; 
    case 4: return { {"VALUE", boost::any(boost::get<QString>(_vmap))} }; 
  }
}
ClCValue::ClCValue_iterator ClCValue::begin() const {
  ClCValue_iterator it;
  it._me = this;
  switch (_vmap.which()) {
    case 0: it._isend = true; it._isfirst = false; break;
    case 1: it._isend = false; it._isfirst = true; it._it = boost::get<const QMap<QString, boost::any>&>(_vmap).begin(); break;
    case 2: it._isend = false; it._isfirst = true; break;
    case 3: it._isend = false; it._isfirst = true; break;
    case 4: it._isend = false; it._isfirst = true; break;
  }
  return it;
}

ClCValue::ClCValue_iterator ClCValue::end() const {
  ClCValue_iterator it;
  it._isfirst = false;
  it._isend = true;
  it._me = this;
  return it;
}

ClCValue::ClCValue_iterator& ClCValue::ClCValue_iterator::operator++() {
  if (_isend) return *this;
  switch (_me->_vmap.which()) {
    case 1: {
      ++_it;
      if (boost::get<const QMap<QString, boost::any>&>(_me->_vmap).end() == _it) {
        _isend = true;
      }
      break;
    }
    default: _isend = true;

  }
  return *this;
}
const QString& ClCValue::ClCValue_iterator::key() {
  static QString null_const;
  static QString value_const = "VALUE";
  if (_isend) return null_const;
  if (_me->_vmap.which() != 1) return value_const;
  return _it.key();
}

boost::any ClCValue::ClCValue_iterator::valueAny() {
  static boost::any fail_dummy;
  if (_isend) {
    EC_THROW("Access end iterator");
    return fail_dummy;
  }
  switch (_me->_vmap.which()) {
    default:
    case 0: return {};
    case 1: return *_it;
    case 2: return boost::get<int>(_me->_vmap);
    case 3: return boost::get<double>(_me->_vmap);
    case 4: return boost::get<QString>(_me->_vmap);
  }
  return {};
}

bool ClCValue::ClCValue_iterator::operator==(const ClCValue::ClCValue_iterator& other) const {
  if (_me != other._me) return false;
  if (_isend && other._isend) return true;
  if (_isfirst && other._isfirst && _me == other._me) return true;
  if (_me->_vmap.which() == 1 && other._me->_vmap.which() == 1 && _it == other._it) return true;
  return false;
}
bool ClCValue::ClCValue_iterator::operator!=(const ClCValue::ClCValue_iterator& other) const {
  return !(*this == other);
}

void ClCValue::shiftDepths(double delta) {
  if (hasValue<ClDepthPointVector>(POINTSET_SCL)) {
    auto dpv = value<ClDepthPointVector>(POINTSET_SCL);
    dpv.shiftDepth(delta);
    setValue<ClDepthPointVector>(POINTSET_SCL, dpv);
  }
  if (hasValue<ClDepthPointVector>(POINTSET_REL)) {
    auto dpv = value<ClDepthPointVector>(POINTSET_REL);
    dpv.shiftDepth(delta);
    setValue<ClDepthPointVector>(POINTSET_REL, dpv);
  }
  if (hasValue<ClDepthPointVector>(POINTSET_ABS)) {
    auto dpv = value<ClDepthPointVector>(POINTSET_ABS);
    dpv.shiftDepth(delta);
    setValue<ClDepthPointVector>(POINTSET_ABS, dpv);
  }
  if (hasValue<ClDepthPointVector>(POINTSET_PHI)) {
    auto dpv = value<ClDepthPointVector>(POINTSET_PHI);
    dpv.shiftDepth(delta);
    setValue<ClDepthPointVector>(POINTSET_PHI, dpv);
  }
}

#if 0
void ClCValue::dumpContents() {
  switch (_vmap.which()) {
  case 0: qDebug() << "VALUE=NULL"; break;
  case 2: qDebug() << "VALUE=(int)" << valueBase<int>(); break;
  case 3: qDebug() << "VALUE=(double)" << valueBase<double>(); break;
  case 4: qDebug() << "VALUE=(QString)" << valueBase<QString>(); break;
  case 1: {
    auto mp = boost::get<QMap<QString, boost::any>>(&_vmap);
    for (auto it = mp->begin(); it != mp->end(); ++it) {
      QString str = it.key() + "=(" + it.value().type().name() + ")";
      if (it.value().type() == typeid(double)) {
        str += QString::number(anyc<double>(it.value()));
      }
      else if (it.value().empty()) { 
        str += "NULL";
      } 
      qDebug() << str;
    }
    break;
  }
  }
}
#endif

#if 0
double ClIntervalMapIterator::topDepth() const {
  if (dmit == _ivmap->_dm->end())  throw std::out_of_range("ClIntervalMapIterator");
  if (_ivmap->_dm->isInterval()) {
    if (boost::next(dmit, 1) == _ivmap->_dm->end() ) throw std::out_of_range("ClIntervalMapIterator");
  }
  return dmit->first;
}

double ClIntervalMapIterator::bottomDepth() const {
  if (dmit == _ivmap->_dm->end())  throw std::out_of_range("ClIntervalMapIterator");
  if (_ivmap->_dm->isInterval()) {
    if (boost::next(dmit, 1) == _ivmap->_dm->end()) throw std::out_of_range("ClIntervalMapIterator");
    return boost::next(dmit, 1)->first;
  }
  else if (_ivmap->_dm->dataType() == CDT_Image) {
    if (dmit->second.type()==typeid(ClColumnImage)) {
      ClColumnImage ci = anyc<ClColumnImage>(dmit->second);
      return ci._bottomdepth;
    }
  }
  else if (_ivmap->_dm->dataType() == CDT_Notes) {
    if (dmit->second.type() == typeid(ClNotesEntry)) {
      ClNotesEntry ne = anyc<ClNotesEntry>(dmit->second);
      return ne.bottomdepth;
    }

  }
  return dmit->first;
}

boost::any ClIntervalMapIterator::value() const {
  if (dmit == _ivmap->_dm->end() )  throw std::out_of_range("ClIntervalMapIterator");
  if (_ivmap->_dm->isInterval()) {
    if (boost::next(dmit, 1) == _ivmap->_dm->end()) throw std::out_of_range("ClIntervalMapIterator");
    return dmit->second;
  }
  return dmit->second;
}

boost::any ClIntervalMapIterator::namedProperty(const QString &) const {
  return value();
}

ClIntervalMapIterator& ClIntervalMapIterator::operator++() {
  dmit++;
  if (_ivmap->_dm->isInterval()) {
    while (dmit != _ivmap->_dm->end() && dmit->second.empty()) ++dmit;
  }
  return *this;
}

ClIntervalMapIterator& ClIntervalMapIterator::operator--() {
  --dmit;
  if (_ivmap->_dm->isInterval()) {
    while (!dmit->second.empty()) --dmit;
  }
  return *this;
}

bool ClIntervalMapIterator::operator==(const ClIntervalMapIterator& other) const {
  return _ivmap == other._ivmap && dmit == other.dmit;
}

bool ClIntervalMapIterator::operator!=(const ClIntervalMapIterator& other) const {
  return !(*this==other);
}

ClIntervalMapIterator ClIntervalMap::begin() const {
  auto dmit = _dm->begin();
  while (dmit != _dm->end() && dmit->second.empty()) ++dmit;
  return ClIntervalMapIterator(this, dmit);
}

ClIntervalMapIterator ClIntervalMap::end() const {
  return ClIntervalMapIterator(this, _dm->end());
}

#endif
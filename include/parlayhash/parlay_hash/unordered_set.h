#ifndef PARLAY_UNORDERED_SET_
#define PARLAY_UNORDERED_SET_

#include <functional>
#include <optional>
#include "parlay_hash.h"
#include <utils/epoch.h>

namespace parlay {

  // entries just contain a key
  template <typename K_, class Hash_ = std::hash<K_>, class KeyEqual_ = std::equal_to<K_>>
  struct SetData {
    using K = K_;
    using Hash = Hash_;
    using KeyEqual = KeyEqual_;
    using value_type = K;
    static const K& get_key(const value_type& x) { return x;}
  };

  // Generic unordered_set that can be used with direct or indirect
  // entries depending on the template argument.
  template <typename Entries>
  struct unordered_set_internal {
    using set = parlay_hash<Entries>;

    Entries entries_;
    set m;

    using Entry = typename Entries::Entry;
    using K = typename Entries::DataS::K;
    using key_type = K;
    using value_type = K;
    using iterator = typename set::Iterator;

    static constexpr auto true_f = [] (const Entry& kv) {return true;};
    static constexpr auto identity = [] (const Entry& kv) {return kv;};

    unordered_set_internal(long n, bool clear_at_end = default_clear_at_end)
      : entries_(Entries(clear_at_end)),
	m(set(n, &entries_, clear_at_end)) {}
    
    iterator begin() { return m.begin();}
    iterator end() { return m.end();}
    bool empty() { return size() == 0;}
    bool max_size() { return (1ul << 47)/sizeof(Entry);}
    void clear() { m.clear_buckets();}
    long size() { return m.size();}

    template <typename F = decltype(identity)>
    auto entries(const F& f = identity) { return m.entries(f);}
    long count(const K& k) { return (contains(k)) ? 1 : 0; }
    bool contains(const K& k) { return find(k, true_f).has_value();}

    bool Find(const K& k) { return m.Find(Entry::make_key(k), true_f).has_value(); }
    bool Insert(const K& key) 
    {
      auto k = Entry::make_key(key);
      return !m.Insert(k, [&] {return entries_.make_entry(k, key);}, true_f).has_value();
    }

    bool Remove(const K& k)
    { return m.Remove(Entry::make_key(k), true_f).has_value(); }

    iterator find(const K& k) { return m.find(k); }

    std::pair<iterator,bool> insert(const value_type& entry) {
      return m.insert(entries_.make_entry(make_key(entry.first), entry)); }

    iterator erase(iterator pos) { return m.erase(pos); }
    size_t erase(const K& k) { return m.erase(k); }

  };

  template <typename K, class Hash = std::hash<K>, class KeyEqual = std::equal_to<K>>
  using parlay_unordered_set_direct = unordered_set_internal<DirectEntries<SetData<K, Hash, KeyEqual>>>;

  template <typename K, class Hash = std::hash<K>, class KeyEqual = std::equal_to<K>>
  using parlay_unordered_set_indirect = unordered_set_internal<IndirectEntries<SetData<K, Hash, KeyEqual>>>;

  template <typename K, class Hash = std::hash<K>, class KeyEqual = std::equal_to<K>>
  using parlay_unordered_set = std::conditional_t<std::is_trivially_copyable_v<K>,
						  parlay_unordered_set_direct<K, Hash, KeyEqual>,
						  parlay_unordered_set_indirect<K, Hash, KeyEqual>>;
}  // namespace parlay
#endif  // PARLAY_BIGATOMIC_HASH_LIST


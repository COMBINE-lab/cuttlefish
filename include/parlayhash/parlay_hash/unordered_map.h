// Initial Author: Guy Blelloch
// Developed as part of the flock library
// 
// A growable unordered_map using a hash table designed for scalability to large number of threads, and
// for high contention.  On a key type K and value type V it supports:
//
//   unordered_map<K, V, Hash=std::hash<K>, Equal=std::equal_to<K>>(n) :
//   constructor for table of initial size n
//
//   Find(const K&) -> std::optional<V> :
//   returns value if key is found, and otherwise returns nullopt
//
//   Insert(const K&, const V&) -> std::optional<V> :
//   if key not in the table it inserts the key with the given value
//   and returns nullopt, otherwise it does not modify the table and
//   returns the old value.
//
//   Remove(const K&) -> std::optional<V> :
//   if key is in the table it removes the entry and returns its value.
//   otherwise it does nothing and returns nullopt.
//
//   size() -> long : returns the size of the table.  Not linearizable with
//   the other functions, and takes time proportional to the table size.
//  
//   clear() -> void : clears the table so its size is 0.
//
//   for_each(F f) : applies functor f to each entry of the table.
//   f should be of type (const std::pair<K,V>&) -> void

#ifndef PARLAY_UNORDERED_MAP_
#define PARLAY_UNORDERED_MAP_

#include <functional>
#include <optional>
#include "parlay_hash.h"

namespace parlay {

  // entries contain a key
  template <typename K_, typename V_, class Hash_ = std::hash<K_>, class KeyEqual_ = std::equal_to<K_>>
  struct MapData {
    using K = K_;
    using V = V_;
    using Hash = Hash_;
    using KeyEqual = KeyEqual_;
    using value_type = std::pair<K,V>;
    static const K& get_key(const value_type& x) { return x.first;}
  };

  // Generic unordered_map that can be used with direct or indirect
  // entries depending on the template argument.
  template <typename Entries>
  struct unordered_map_internal {
    using map = parlay_hash<Entries>;

    Entries entries_;
    map m;

    using Entry = typename Entries::Entry;
    using K = typename Entries::DataS::K;
    using V = typename Entries::DataS::V;
    using key_type = K;
    using mapped_type = V;
    using value_type = std::pair<K, V>;
    using iterator = typename map::Iterator;

    static constexpr auto true_f = [] (const Entry& kv) {return true;};
    static constexpr auto identity = [] (const Entry& kv) {return kv;};
    static constexpr auto get_value = [] (const value_type& kv) {return kv.second;};

    unordered_map_internal(long n, bool clear_at_end = default_clear_at_end)
      : entries_(Entries(clear_at_end)),
	m(map(n, &entries_, clear_at_end)) {}
    
    iterator begin() { return m.begin();}
    iterator end() { return m.end();}
    bool empty() { return size() == 0;}
    bool max_size() { return (1ul << 47)/sizeof(Entry);}
    void clear() { m.clear_buckets();}
    long size() { return m.size();}

    template <typename F = decltype(identity)>
    //auto entries(const F& f = identity) { return m.entries(f);}
    long count(const K& k) { return (contains(k)) ? 1 : 0; }
    bool contains(const K& k) { return find(k, true_f).has_value();}

    template <typename F = decltype(get_value)>
    auto Find(const K& k, const F& f = get_value)
      -> std::optional<typename std::result_of<F(value_type)>::type>
    {
      auto g = [&] (const Entry& e) {return f(e.get_entry());};
      return m.Find(Entry::make_key(k), g);
    }

    auto Insert(const K& key, const V& value) -> std::optional<mapped_type>
    {
      auto k = Entry::make_key(key);
      auto g = [&] (const Entry& e) {return get_value(e.get_entry());};
      return m.Insert(k, [&] {return entries_.make_entry(k, value_type(key, value));}, g);
    }

    template <typename F>
    auto Upsert(const K& key, const F& f) -> std::optional<mapped_type>
    {
      auto k = Entry::make_key(key);
      auto g = [&] (const Entry& e) {return get_value(e.get_entry());};
      auto constr = [&] (const std::optional<Entry>& e) -> Entry {
		      if (e.has_value())
			return entries_.make_entry(k, value_type(key, f(std::optional(get_value((*e).get_entry())))));
		      return entries_.make_entry(k, value_type(key, f(std::optional<V>())));
		    };
      return m.Upsert(k, constr, g);
    }

    template <typename F>
    auto Insert(const K& key, const V& value, const F& f)
      -> std::optional<typename std::result_of<F(value_type)>::type>
    {
      auto k = Entry::make_key(key);
      auto g = [&] (const Entry& e) {return f(e.get_entry());};
      return m.Insert(k, [&] {return entries_.make_entry(k, value_type(key, value));}, g);
    }

    auto Remove(const K& k) -> std::optional<mapped_type>
    {
      auto g = [&] (const Entry& e) {return get_value(e.get_entry());};
      return m.Remove(Entry::make_key(k), g);
    }

    template <typename F>
    auto Remove(const K& k, const F& f)
      -> std::optional<typename std::result_of<F(value_type)>::type>
    {
      auto g = [&] (const Entry& e) {return f(e.get_entry());};
      return m.Remove(Entry::make_key(k), g);
    }

    iterator find(const K& k) { return m.find(k); }

    std::pair<iterator,bool> insert(const value_type& entry) {
      auto k = Entry::make_key(entry.first);
      return m.insert(k, [&] {return entries_.make_entry(k, entry);});}

    iterator erase(iterator pos) { return m.erase(pos); }
    size_t erase(const K& k) { return m.erase(k); }

  };

  // Entries are stored directly in the bucket, avoiding a cache miss
  // for indirection.  Entries can be moved by updates even on
  // different keys.
  template <typename K, typename V, class Hash = std::hash<K>, class KeyEqual = std::equal_to<K>>
  using parlay_unordered_map_direct = unordered_map_internal<DirectEntries<MapData<K, V, Hash, KeyEqual>>>;

  // Entries are stored indirectly through a pointer.  Pointers to
  // entries wil remain valid until the entry is upserted or deleted
  // (an upsert can be though of as a deletion followed by an
  // insersion).
  template <typename K, typename V, class Hash = std::hash<K>, class KeyEqual = std::equal_to<K>>
  using parlay_unordered_map_indirect = unordered_map_internal<IndirectEntries<MapData<K, V, Hash, KeyEqual>>>;

  // specialization of unordered_map to use either direct or indirect
  // entries depending on whether K and V are trivially copyable.
  template <typename K, typename V, class Hash = std::hash<K>, class KeyEqual = std::equal_to<K>>
  using parlay_unordered_map = std::conditional_t<std::is_trivially_copyable_v<K> &&
						  std::is_trivially_copyable_v<V>,
						  parlay_unordered_map_direct<K,V,Hash,KeyEqual>,
						  parlay_unordered_map_indirect<K,V,Hash,KeyEqual>>;
}  // namespace parlay
#endif  // PARLAY_BIGATOMIC_HASH_LIST

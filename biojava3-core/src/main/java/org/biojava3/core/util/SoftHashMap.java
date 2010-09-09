/*
 *                    PDB web development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 *
 * Created on Mar 26, 2009
 * Created by ap3
 *
 */

package org.biojava3.core.util;

import java.lang.ref.ReferenceQueue;
import java.lang.ref.SoftReference;
import java.util.AbstractMap;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;


/** A in memory cache using soft references. (can be garbage collected)
 * 
 * This code is based on: http://java-interview-faqs.blogspot.com/2008/09/building-faster-and-efficient-cache.html 
 * */


public class SoftHashMap<K, V> extends AbstractMap<K, V> {

   public static final boolean DEBUG = false;


   public static final int DEFAULT_LIMIT = 1;

   /** The internal HashMap that stores SoftReference to actual data. */

   private final Map<K, SoftReference<V>> map = new HashMap<K, SoftReference<V>>();

   /** Maximum Number of references you dont want GC to collect. */

   private final int max_limit;

   /** The FIFO list of hard references, order of last access. */

   private final LinkedList<V> hardCache = new LinkedList<V>();

   /** Reference queue for cleared SoftReference objects. */

   private final ReferenceQueue<V> queue = new ReferenceQueue<V>();

   public SoftHashMap() {

      this(1000);

   }

   public SoftHashMap(int hardSize) {

      max_limit = hardSize;

   }

   public V get(Object key) {

      V result = null;

      // We get the SoftReference represented by that key

      SoftReference<V> soft_ref = map.get(key);

      if (soft_ref != null) {

         try {

            // From the SoftReference we get the value, which can be

            // null if it was not in the map, or it was removed in

            // the clearGCCollected() method defined below

            result = soft_ref.get();

            if (result == null) {

               // If the value has been garbage collected, remove the

               // entry from the HashMap.

               map.remove(key);

            } else {

               // We now add this object to the beginning of the hard

               // reference queue. One reference can occur more than

               // once, because lookups of the FIFO queue are slow, so

               // we don't want to search through it each time to remove

               // duplicates.

               synchronized (hardCache){
                  hardCache.addFirst(result);

                  if (hardCache.size() > max_limit) {

                     // Remove the last entry if list greater than MAX_LIMIT

                     hardCache.removeLast();

                  }
               }

            }
         } catch (Exception e){
            e.printStackTrace();
         }

      }

      return result;

   }



   /**

    * We define our own subclass of SoftReference which contains not only the

    * value but also the key to make it easier to find the entry in the HashMap

    * after it's been garbage collected.

    */

   private static class SoftValue<K, V> extends SoftReference<V> {

      private final Object key; // always make data member final



      /**

       * Did you know that an outer class can access private data members and

       * methods of an inner class? I didn't know that! I thought it was only

       * the inner class who could access the outer class's private

       * information. An outer class can also access private members of an

       * inner class inside its inner class.

       */

      private SoftValue(V k, K key, ReferenceQueue<? super V> q) {

         super(k, q);

         this.key = key;

      }

   }



   /**

    * Here we go through the ReferenceQueue and remove garbage collected

    * SoftValue objects from the HashMap by looking them up using the

    * SoftValue.key data member.

    */

   @SuppressWarnings("unchecked") // every Reference in queue is stored as a SoftValue
   private void clearGCCollected() {

      SoftValue<K, V> sv;

      while ((sv = (SoftValue<K, V>) queue.poll()) != null) {

         map.remove(sv.key); // we can access private data!

      }

   }



   /**

    * Here we put the key, value pair into the HashMap using a SoftValue

    * object.

    */

   public synchronized V put(K key, V value) {

      clearGCCollected();

      if ( DEBUG)
         System.out.println("putting " + key + " on cache. size: " + size());
      
      map.put(key, new SoftValue<K, V>(value, key, queue));
      
      return value;

   }



   public V remove(Object key) {

      clearGCCollected();
      if ( DEBUG)
          System.out.println("removing " + key + " from cache. size: " + size());
      return map.remove(key).get();

   }



   public void clear() {

      synchronized (hardCache){
         hardCache.clear();
      }

      clearGCCollected();
      if ( DEBUG)
          System.out.println("clearing cache");
      map.clear();

   }



   public int size() {

      clearGCCollected();

      return map.size();

   }



   public Set<Map.Entry<K, V>> entrySet() {

      throw new UnsupportedOperationException();

   }

}
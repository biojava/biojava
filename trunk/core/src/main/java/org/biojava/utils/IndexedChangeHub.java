/*
 *                    BioJava development code
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
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */


package org.biojava.utils;

import java.lang.ref.Reference;
import java.lang.ref.ReferenceQueue;
import java.lang.ref.WeakReference;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * implements Changeable support with a ChangeHub that
 * stores ChangeListener by key.
 *
 * @author Thomas Down (original implementation)
 * @author David Huen (refactoring)
 * @since 1.3
 */
public abstract class IndexedChangeHub implements ChangeHub
{
    private ReferenceQueue refQueue;
    private Map listeners;

    public IndexedChangeHub()
    {
        refQueue = new ReferenceQueue();
        listeners = new HashMap();
    }

    // queue cleanup
    // the references are now safe

    abstract protected boolean isMyChangeEvent(ChangeEvent cev, IndexedChangeHub.ListenerMemento lm);

    public void addListener(Object key, ChangeListener listener, ChangeType ct)
    {
        diddleQueue();
        List listenerList = (List) listeners.get(key);
        if (listenerList == null) {
            listenerList = new ArrayList();
            listeners.put(key, listenerList);
        }
        listenerList.add(new ListenerMemento(ct, new ListenerReference(key, listener, refQueue)));
    }

    public void removeListener(Object key, ChangeListener listener, ChangeType ct)
    {
        List listenerList = (List) listeners.get(key);
        if (listenerList != null) {
            for (Iterator i = listenerList.iterator(); i.hasNext(); ) {
                ListenerMemento lm = (ListenerMemento) i.next();
                if (ct == lm.type && listener.equals(lm.listener.get())) {
                    lm.listener.clear();
                    i.remove();
                    return;
                }
            }
        }
    }

    public void firePreChange(Object key, ChangeEvent cev)
        throws ChangeVetoException
    {
        List listenerList = (List) listeners.get(key);
        if (listenerList != null) {
            for (Iterator i = listenerList.iterator(); i.hasNext(); ) {
                ListenerMemento lm = (ListenerMemento) i.next();
                if (isMyChangeEvent(cev, lm)) {
                    ChangeListener cl = (ChangeListener) lm.listener.get();
                    if (cl != null) {
                        cl.preChange(cev);
                    }
                }
            }
        }

        // in the original version, there was the possibility of firing a
        // ChangeEvent for the parent class here.  This is not feasible
        // in this implementation.  The child must override this method
        // and fire it themselves.

    }

    public void firePostChange(Object key, ChangeEvent cev)
    {
        List listenerList = (List) listeners.get(key);
        if (listenerList != null) {
            for (Iterator i = listenerList.iterator(); i.hasNext(); ) {
                ListenerMemento lm = (ListenerMemento) i.next();
                if (isMyChangeEvent(cev, lm)) {
                    ChangeListener cl = (ChangeListener) lm.listener.get();
                    if (cl != null) {
                        cl.postChange(cev);
                    }
                }
            }
        }
        // in the original version, there was the possibility of firing a
        // ChangeEvent for the parent class here.  This is not feasible
        // in this implementation.  The child must override this method
        // and fire it themselves.

    }

    protected void diddleQueue()
    {
        Reference ref;
        while ((ref = refQueue.poll()) != null) {
            List listenerList = (List) listeners.get(((ListenerReference) ref).getKey());
            if (listenerList != null) {
                for (Iterator i = listenerList.iterator(); i.hasNext(); ) {
                    ListenerMemento lm = (ListenerMemento) i.next();
                    if (lm.listener == ref) {
                        i.remove();
                        break;
                    }
                }
            }
        }
    }

    private class ListenerReference extends WeakReference {
        private Object key;

        public ListenerReference(Object key, Object ref) {
            super(ref);
            this.key = key;
        }

        public ListenerReference(Object key, Object ref, ReferenceQueue queue) {
            super(ref, queue);
            this.key = key;
        }

        public Object getKey() {
            return key;
        }
    }

    protected class ListenerMemento {
        public final ChangeType type;
        public final Reference listener;

        public ListenerMemento(ChangeType type, Reference listener) {
            this.type = type;
            this.listener = listener;
        }
    }

}


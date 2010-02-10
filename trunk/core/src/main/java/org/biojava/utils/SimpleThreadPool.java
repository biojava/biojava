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

import java.util.LinkedList;

/**
 * <p><code>SimpleThreadPool</code> is a basic implementation of
 * <code>ThreadPool</code> for use where we don't wish to introduce a
 * dependency on a 3rd-party pool. In general, objects which require a
 * pool should only use the interface and parameterize such that other
 * implementations may be dropped in in place of this one, possibly
 * using this one as a fallback.</p>
 *
 * <p>This class offers a service for running <code>Runnable</code>s
 * using multiple threads, the number of which is specified in the
 * constructor. <code>Runnable</code>s are queued in a simple FIFO
 * queue. The worker threads wait on the queue when it is empty and
 * are notified when a new <code>Runnable</code> is submitted.</p>
 *
 * <p>This implementation will prevent an application from exiting
 * until <code>stopThreads()</code> is called unless the pool contains
 * daemon threads.</p>
 *
 * @author Keith James
 * @since 1.3
 */
public class SimpleThreadPool implements ThreadPool
{
    protected PooledThread [] threads;
    protected int priority;

    private LinkedList queue;
    private boolean daemon;
    private boolean waiting;
    private boolean stopped;

    /**
     * Creates a new <code>SimpleThreadPool</code> containing 4
     * non-daemon threads and starts them. The threads have priority
     * Thread.NORM_PRIORITY.  Because threads are non-deamon you will need
     * to call stopThreads() to terminate them.
     */
    public SimpleThreadPool()
    {
        this(4, false);
    }

    /**
     * Creates a new <code>SimpleThreadPool</code> containing the
     * specified number of threads and starts them. The threads have
     * priority Thread.NORM_PRIORITY.
     *
     * @param threadCount an <code>int</code> thread count.
     * @param daemon a <code>boolean</code> indicating whether the
     * threads should be daemons. If threads are non-deamon you will need
     * to call stopThreads() to terminate them.
     */
    public SimpleThreadPool(int threadCount, boolean daemon)
    {
        this(threadCount, daemon, Thread.NORM_PRIORITY);
    }

    /**
     * Creates a new <code>SimpleThreadPool</code> containing the
     * specified number of threads and starts them.
     *
     * @param threadCount an <code>int</code> thread count.
     * @param daemon a <code>boolean</code> indicating whether the
     * threads should be daemons. If threads are non-deamon you will need
     * to call stopThreads() to terminate them.
     * @param priority an <code>int</code> priority for the threads.
     */
    public SimpleThreadPool(int threadCount, boolean daemon, int priority)
    {
        this.daemon = daemon;
        this.priority = priority;
        queue = new LinkedList();
        threads = new PooledThread[threadCount];
        stopped = true;
        waiting = false;
        startThreads();
    }

    public void addRequest(Runnable task)
    {
        if (waiting || stopped)
            throw new IllegalStateException("Thread pool has been closed to new requests");

        synchronized(queue)
        {
            queue.add(task);
            // Notify threads blocked in nextRequest()
            queue.notifyAll();
        }
    }

    public void startThreads()
    {
        if (! stopped)
            throw new IllegalStateException("Thread pool is already started");

        stopped = false;

        synchronized(threads)
        {
            for (int i = 0; i < threads.length; i++)
            {
                threads[i] = new PooledThread();
                if (daemon)
                    threads[i].setDaemon(true);
                threads[i].setPriority(priority);
                threads[i].start();
            }
        }
    }

    /**
     * Waits for all working threads to return and then stops them. If the
     * thread pool contains non-daemon threads you will have to call this method
     * to make your program return.
     * @throws IllegalStateException if the pool is already stopped.
     */
    public void stopThreads()
    {
        if (stopped)
            throw new IllegalStateException("Thread pool has already been stopped");

        stopped = true;

        synchronized(queue)
        {
            // Ensure working threads return and die
            while (threadsAlive() > 0)
            {
                try
                {
                    queue.wait(500);
                    queue.notifyAll();
                }
                catch (InterruptedException ie) { }
            }
        }
    }

    public void waitForThreads()
    {
        if (stopped)
            throw new IllegalStateException("Thread pool has been stopped");

        waiting = true;

        synchronized(threads)
        {
            // Ensure queue gets emptied and all work is done
            while (! queue.isEmpty() || threadsWorking() > 0)
            {
                try
                {
                    threads.wait();
                }
                catch (InterruptedException ie) { }
            }
        }

        waiting = false;
    }

    /**
     * <code>threadsWorking</code> returns the number of threads
     * currently performing work.
     *
     * @return an <code>int</code>.
     */
    public int threadsWorking()
    {
        int workingCount = 0;

        synchronized(threads)
        {
            for (int i = 0; i < threads.length; i++)
                if (threads[i].working)
                    workingCount++;
        }

        return workingCount;
    }

    /**
     * <code>threadsIdle</code> returns the number of threads
     * currently waiting for work.
     *
     * @return an <code>int</code>.
     */
    public int threadsIdle()
    {
        return threads.length - threadsWorking();
    }

    /**
     * <code>requestsQueued</code> returns the number of
     * <code>Runnable</code>s currently queued.
     *
     * @return an <code>int</code>.
     */
    public int requestsQueued()
    {
        return queue.size();
    }

    /**
     * <code>threadsAlive</code> returns the number of threads
     * currently alive.
     *
     * @return an <code>int</code>.
     */
    protected int threadsAlive()
    {
        int aliveCount = 0;

        synchronized(threads)
        {
            for (int i = 0; i < threads.length; i++)
                if (threads[i].isAlive())
                    aliveCount++;
        }

        return aliveCount;
    }

    /**
     * <code>nextRequest</code> gets the next <code>Runnable</code>
     * from the queue. This method blocks if the queue is empty and
     * the pool has not stopped. If the pool has stopped it returns
     * null.
     *
     * @return a <code>Runnable</code> or null if the pool has been
     * stopped.
     */
    protected Runnable nextRequest()
    {
        synchronized(queue)
        {
            try
            {
                while (! stopped && queue.isEmpty())
                    queue.wait();
            }
            catch (InterruptedException ie) { }

            if (stopped)
                return null;
            else
                return (Runnable) queue.removeFirst();
        }
    }

    /**
     * <code>PooledThread</code> is a thread class which works within
     * the pool. It sets its boolean flag true when working,
     * synchronizing this on the array which contains all the
     * <code>PooledThread</code>s.
     */
    private class PooledThread extends Thread
    {
        boolean working = false;

        public void run()
        {
            while (true)
            {
                Runnable task = nextRequest();

                // If the pool is stopped the queue returns null and
                // the thread dies
                if (task == null)
                    break;

                // Synchronize on thread array to update state
                synchronized(threads)
                {
                    working = true;
                }

                task.run();

                synchronized(threads)
                {
                    working = false;
                    threads.notify();
                }
            }
        }
    }
}

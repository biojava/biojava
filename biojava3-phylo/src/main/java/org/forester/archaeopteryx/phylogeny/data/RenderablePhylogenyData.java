// $Id: RenderablePhylogenyData.java,v 1.3 2009/12/30 04:33:45 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx.phylogeny.data;

import java.awt.Dimension;
import java.awt.Graphics2D;

import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.data.PhylogenyData;

public interface RenderablePhylogenyData extends PhylogenyData {

    public Dimension getOriginalSize();

    public Object getParameter();

    public Dimension getRenderingSize();

    /**
     * This can be used to render phylogeny data as graphics (for example,
     * display of the domain structure). In most Renderable implementations this
     * will do nothing (i.e. just return).
     * 
     * @param g
     *            the Graphics to render to
     */
    public void render( final double x, final double y, final Graphics2D g, final TreePanel tree_panel, boolean to_pdf );

    public void setParameter( final double parameter );

    public void setRenderingHeight( final double rendering_height );
}

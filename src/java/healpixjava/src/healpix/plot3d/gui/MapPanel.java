/*
 * HEALPix Java code supported by the Gaia project.
 * Copyright (C) 2006-2011 Gaia Data Processing and Analysis Consortium
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */
package healpix.plot3d.gui;

import healpix.core.dm.HealpixMap;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextArea;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Panel to Display Map Attributes
 * 
 * @author ejoliet
 * @version $Id: MapPanel.java,v 1.1 2008/04/25 14:44:51 healpix Exp $
 * 
 */
public class MapPanel extends Panel implements ActionListener {
	private static final long serialVersionUID = 1L;

	protected HealpixMap map;

	Label id;

	Label nside;

	Label fwhm;

	Label creator;

	Label version;

	Label extname;

	Label status;

	TextArea notes;

	/**
	 * initialiser
	 */
	protected void init() {
		id = new Label(" ", Label.RIGHT);
		id.setForeground(Color.black);
		nside = new Label("    ", Label.RIGHT);
		nside.setForeground(Color.black);
		fwhm = new Label("    ", Label.RIGHT);
		fwhm.setForeground(Color.black);
		creator = new Label("                    ", Label.RIGHT);
		creator.setForeground(Color.black);
		version = new Label("         ", Label.RIGHT);
		version.setForeground(Color.black);
		extname = new Label("         ", Label.RIGHT);
		extname.setForeground(Color.black);
		notes = new TextArea("", 5, 25, TextArea.SCROLLBARS_BOTH);
		notes.setForeground(Color.black);
		Label idl = new Label("ID:");
		idl.setForeground(Color.blue);
		Label nsidel = new Label("NSIDE:");
		nsidel.setForeground(Color.blue);
		Label fwhml = new Label("FWHM:");
		fwhml.setForeground(Color.blue);
		Label creatorl = new Label("Creator:");
		creatorl.setForeground(Color.blue);
		Label verl = new Label("Version:");
		verl.setForeground(Color.blue);
		Label extl = new Label("Extname:");
		extl.setForeground(Color.blue);
		Label classifl = new Label("Class:");
		classifl.setForeground(Color.blue);
		setLayout(new BorderLayout(0, 0));
		Label title = new Label("Map Details", Label.CENTER);
		title.setForeground(Color.blue);
		add("North", title);
		add("South", notes);
		Panel data = new Panel();
		data.setLayout(new GridLayout(6, 2, 0, 0));
		data.add(idl);
		data.add(id);
		data.add(nsidel);
		data.add(nside);
		data.add(fwhml);
		data.add(fwhm);
		data.add(creatorl);
		data.add(creator);
		data.add(verl);
		data.add(version);
		data.add(extl);
		data.add(extname);
		add("Center", data);
		setBackground(Color.white);
	}

	/**
	 * Set Status output label
	 */
	public void setStatus(Label status) {
		this.status = status;
	}

	/**
	 * Set Map
	 */
	public void setMap(HealpixMap map) {
		this.map = map;
		refresh();
	}

	protected void refresh() {
		// System.out.println("Refresh map panel with map "+map.getId());
		if (map != null) {
			// id.setText(map.getId());
			// fwhm.setText(""+map.getFWHM());
			nside.setText("" + map.nside());
			// creator.setText(map.getCreator());
			// version.setText(map.getVersion());
			// extname.setText(map.getExtname());
			// notes.setText(map.getNotes());
		}
	}

	public MapPanel() {
		init();
	}

	public HealpixMap getMap() {
		return map;
	}

	/**
	 * This is required to implement ActionListener - this gets called when
	 * someone hits a button
	 */
	public void actionPerformed(ActionEvent e) {
	}
}

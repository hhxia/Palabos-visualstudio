// LBM.h : main header file for the LBM application
//
#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"       // main symbols


// CLBMApp:
// See LBM.cpp for the implementation of this class
//

class CLBMApp : public CWinApp
{
public:
	CLBMApp();


// Overrides
public:
	virtual BOOL InitInstance();

// Implementation
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CLBMApp theApp;
#pragma once
#if defined(_MSC_VER)

#if defined(_DEBUG)
#include <crtdbg.h>
#endif // defined(_DEBUG)

namespace testing
{
class MSVCMemoryLeakListener : public ::testing::EmptyTestEventListener
{
#if defined(_DEBUG)
	class CRTReportRedirector
	{
	private:
		static int __cdecl crtReportHook(int reportType, char * szMsg, int * retVal)
		{
			if ( retVal )
				*retVal = 0;
	
			std::cerr << szMsg;
	
			if ( retVal )
				*retVal = 0;
	
			return ( reportType == _CRT_ASSERT );
		}
	
	public:
		CRTReportRedirector(void) { InstallHook(); }
		virtual ~CRTReportRedirector(void) { UninstallHook(); }
	
	public:
		static void InstallHook(void) { _CrtSetReportHook2(_CRT_RPTHOOK_INSTALL, &crtReportHook); }
		static void UninstallHook(void) { _CrtSetReportHook2(_CRT_RPTHOOK_REMOVE, &crtReportHook); }
	};
	
public:
	virtual void OnTestStart(const TestInfo &) { _CrtMemCheckpoint(&m_memState); }

	virtual void OnTestEnd(const TestInfo & test_info)
	{
		if ( !test_info.result()->Passed() )
			return;
	
		_CrtMemState stateNow, stateDiff;
		_CrtMemCheckpoint(&stateNow);
	
		if ( !_CrtMemDifference(&stateDiff, &m_memState, &stateNow) )
			return;
	
		{
			CRTReportRedirector _;
			_CrtMemDumpAllObjectsSince(&m_memState);
			_CrtMemDumpStatistics(&stateDiff);
		}
	
		FAIL() << "Memory leak of " << stateDiff.lSizes[1] << " byte(s) detected.";
	}

private:
	_CrtMemState m_memState;
#endif // defined(_DEBUG)
};
}

#endif // defined(_MSC_VER)
